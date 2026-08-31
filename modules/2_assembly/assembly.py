#!/usr/bin/env python3
"""Assemble paired short reads or one long-read metagenomic file."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


DEFAULTS = """Module 2 assembly defaults:
threads=16
memory=51 GB for the assembler (80% of the 64 GB workflow allocation)
min_len=1000
assembler=megahit
megahit_k_min=21
megahit_k_max=141
megahit_k_step=12
megahit_merge_level=20,0.95
metaspades_k_list=21,33,55
long_read_type=required with --long-reads
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
skip_quast=false
keep_temp=false"""

DEFAULT_THREADS = 16
DEFAULT_MEMORY_GB = 51
LONG_READ_TYPES = (
    "pacbio-raw",
    "pacbio-corr",
    "pacbio-hifi",
    "nano-raw",
    "nano-corr",
    "nano-hq",
)


def executable(name: str) -> str:
    result = shutil.which(name)
    if result is None:
        raise RuntimeError(f"Required executable is not on PATH: {name}")
    return result


def run(argv: list[str], *, stdout=None, env: dict[str, str] | None = None) -> None:
    print("[INFO] Executing:", " ".join(argv), flush=True)
    subprocess.run(argv, check=True, stdout=stdout, env=env)


def module_root(project_path: Path) -> Path:
    if (project_path / "2_assembly" / "rm_short_contigs.py").is_file():
        return project_path
    if (project_path / "modules" / "2_assembly" / "rm_short_contigs.py").is_file():
        return project_path / "modules"
    raise FileNotFoundError(f"Cannot locate the flat 2_assembly module from project path: {project_path}")


def validate_kmers(args: argparse.Namespace) -> None:
    if args.assembler == "metaspades":
        for value in args.k_list.split(","):
            try:
                kmer = int(value)
            except ValueError as error:
                raise ValueError(f"metaSPAdes k-mer size is not an integer: {value}") from error
            if kmer % 2 == 0 or kmer >= 128:
                raise ValueError(f"metaSPAdes k-mer size must be odd and below 128: {kmer}")
    if args.assembler == "megahit":
        if args.k_min % 2 == 0 or args.k_max % 2 == 0 or args.k_step % 2 != 0:
            raise ValueError("MEGAHIT k-min/k-max must be odd and k-step must be even")
        if args.k_min < 15 or args.k_max > 255 or args.k_min > args.k_max:
            raise ValueError(f"Invalid MEGAHIT k-mer range: {args.k_min}-{args.k_max}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module assembly --help",
    )
    parser.add_argument("-p", "--project-path", type=Path, required=True, help="METAHICT repository or modules root")
    parser.add_argument("-1", "--r1", type=Path, help="Cleaned short-read shotgun forward FASTQ")
    parser.add_argument("-2", "--r2", type=Path, help="Cleaned short-read shotgun reverse FASTQ")
    parser.add_argument(
        "--long-reads",
        type=Path,
        help="Single long-read shotgun FASTA/FASTQ supplied directly to metaFlye",
    )
    parser.add_argument(
        "--long-read-type",
        choices=LONG_READ_TYPES,
        help="metaFlye input type required with --long-reads",
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Threads passed to the assembler and QUAST")
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=DEFAULT_MEMORY_GB,
        help="MEGAHIT/metaSPAdes memory in GB; metaFlye has no native memory limit",
    )
    parser.add_argument("-l", "--min-len", type=int, default=1000, help="Minimum assembled contig length")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--megahit", dest="assembler", action="store_const", const="megahit", help="Use MEGAHIT")
    group.add_argument("--metaspades", dest="assembler", action="store_const", const="metaspades", help="Use metaSPAdes")
    group.add_argument("--metaflye", dest="assembler", action="store_const", const="metaflye", help="Use metaFlye")
    parser.set_defaults(assembler="megahit")
    parser.add_argument("--k-min", type=int, default=21, help="MEGAHIT minimum k-mer size")
    parser.add_argument("--k-max", type=int, default=141, help="MEGAHIT maximum k-mer size")
    parser.add_argument("--k-step", type=int, default=12, help="MEGAHIT k-mer step size")
    parser.add_argument("--merge-level", default="20,0.95", help="MEGAHIT merge level")
    parser.add_argument("--megahit-extra-args", default="", help="Additional native MEGAHIT arguments")
    parser.add_argument("--k-list", default="21,33,55", help="metaSPAdes comma-separated k-mer list")
    parser.add_argument("--metaspades-extra-args", default="", help="Additional native metaSPAdes arguments")
    parser.add_argument("--metaflye-extra-args", default="", help="Additional native metaFlye arguments")
    parser.add_argument("--tmp-dir", type=Path, help="Temporary directory root")
    parser.add_argument("--skip-quast", action="store_true", help="Skip assembly QC with QUAST")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary assembler files")
    parser.add_argument("--print-defaults", action="store_true", help="Print machine-readable defaults and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.print_defaults:
        print(DEFAULTS)
        return 0
    memory = args.memory
    if min(args.threads, memory, args.min_len) < 1:
        raise ValueError("threads, memory, and minimum contig length must be positive")
    if args.long_reads:
        if args.r1 or args.r2:
            raise ValueError("--long-reads cannot be combined with -1/--r1 or -2/--r2")
        if not args.long_read_type:
            raise ValueError("--long-read-type is required with --long-reads")
        if args.assembler != "metaflye":
            raise ValueError("--long-reads requires --metaflye")
        input_reads = (args.long_reads,)
    else:
        if not args.r1 or not args.r2:
            raise ValueError("paired short-read assembly requires both -1/--r1 and -2/--r2")
        if args.long_read_type:
            raise ValueError("--long-read-type requires --long-reads")
        if args.assembler == "metaflye":
            raise ValueError("--metaflye requires --long-reads and --long-read-type")
        input_reads = (args.r1, args.r2)
    for read in input_reads:
        if not read.is_file():
            raise FileNotFoundError(f"Input read file not found: {read}")
    validate_kmers(args)

    root = module_root(args.project_path)
    short_filter = root / "2_assembly" / "rm_short_contigs.py"
    assembler = {"megahit": "megahit", "metaspades": "metaspades.py", "metaflye": "flye"}[args.assembler]
    assembler_cmd = executable(assembler)
    samtools = executable("samtools")
    quast = None if args.skip_quast else executable("quast.py")
    tmp_root = args.tmp_dir or Path(os.environ.get("METAHICT_TMP_ROOT") or os.environ.get("TMPDIR") or tempfile.gettempdir())
    if not tmp_root.is_dir() or not os.access(tmp_root, os.W_OK):
        raise RuntimeError(f"Temporary directory is not writable: {tmp_root}")

    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    matplotlib = out / ".matplotlib"
    matplotlib.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment["TMPDIR"] = str(tmp_root)
    environment["MPLCONFIGDIR"] = str(matplotlib)
    run_tmp = Path(tempfile.mkdtemp(prefix=f"metahict_{args.assembler}.", dir=tmp_root))
    final = out / "final_assembly.fasta"
    try:
        if args.assembler == "metaspades":
            result_dir = out / "metaspades"
            run([assembler_cmd, "--tmp-dir", str(run_tmp), "-k", args.k_list, "-t", str(args.threads),
                 "-m", str(memory), *shlex.split(args.metaspades_extra_args), "-o", str(result_dir),
                 "-1", str(args.r1), "-2", str(args.r2)], env=environment)
            source = result_dir / "scaffolds.fasta"
            if not source.is_file():
                raise RuntimeError("metaSPAdes did not create scaffolds.fasta")
            with final.open("w") as output_handle:
                run([sys.executable, str(short_filter), str(args.min_len), str(source)], stdout=output_handle, env=environment)
        elif args.assembler == "megahit":
            result_dir = out / "megahit"
            run([assembler_cmd, "-1", str(args.r1), "-2", str(args.r2), "-o", str(result_dir),
                 "--min-contig-len", str(args.min_len), "--k-min", str(args.k_min), "--k-max", str(args.k_max),
                 "--k-step", str(args.k_step), "--merge-level", args.merge_level, "-t", str(args.threads),
                 "-m", str(memory * 1024**3), "--tmp-dir", str(run_tmp),
                 *shlex.split(args.megahit_extra_args)], env=environment)
            source = result_dir / "final.contigs.fa"
            if not source.is_file():
                raise RuntimeError("MEGAHIT did not create final.contigs.fa")
            shutil.copy2(source, final)
        else:
            result_dir = out / "metaflye"
            run([assembler_cmd, f"--{args.long_read_type}", str(args.long_reads), "--meta", "-t", str(args.threads),
                 *shlex.split(args.metaflye_extra_args), "--out-dir", str(result_dir)], env=environment)
            source = result_dir / "assembly.fasta"
            if not source.is_file():
                raise RuntimeError("Flye did not create assembly.fasta")
            shutil.copy2(source, final)
        if not final.is_file() or final.stat().st_size == 0:
            raise RuntimeError("Final assembly is empty")
        if quast:
            run([quast, "-t", str(args.threads), "-o", str(out / "QUAST_out"), str(final)], env=environment)
        run([samtools, "faidx", str(final)], env=environment)
    finally:
        if not args.keep_temp:
            shutil.rmtree(run_tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
