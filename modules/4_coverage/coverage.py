#!/usr/bin/env python3
"""Map paired short reads or single-file long reads and calculate coverage."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


DEFAULT_THREADS = 16
DEFAULT_MEMORY = "25g"


def executable(name: str) -> str:
    result = shutil.which(name)
    if result is None:
        raise RuntimeError(f"Required executable is not on PATH: {name}")
    return result


def memory_default() -> str:
    """Return BBMap heap matching 80% of the 32 GB workflow allocation."""
    return DEFAULT_MEMORY


def run_logged(argv: list[str], log: Path, *, env: dict[str, str] | None = None) -> None:
    print("[INFO] Executing:", " ".join(argv), flush=True)
    with log.open("w") as handle:
        subprocess.run(argv, check=True, stdout=handle, stderr=subprocess.STDOUT, env=env)


def clear_output(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    preserved = {"module.log"}
    for child in directory.iterdir():
        if child.name in preserved:
            continue
        if child.is_dir() and not child.is_symlink():
            shutil.rmtree(child)
        else:
            child.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module coverage --help",
    )
    parser.add_argument("-p", "--project-path", required=True, help="METAHICT repository root")
    parser.add_argument("-1", type=Path, dest="reads1", help="Cleaned short-read shotgun forward FASTQ")
    parser.add_argument("-2", type=Path, dest="reads2", help="Cleaned short-read shotgun reverse FASTQ")
    parser.add_argument(
        "--long-reads",
        type=Path,
        help="Single long-read shotgun FASTA/FASTQ mapped as unpaired reads",
    )
    parser.add_argument("-r", "--reference", type=Path, required=True, help="Assembly FASTA")
    parser.add_argument("-o", "--output", "--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Threads passed to BBMap and SAMtools")
    parser.add_argument(
        "-m",
        "--memory",
        default=DEFAULT_MEMORY,
        help="BBMap Java heap, such as 25g or 25600m",
    )
    parser.add_argument(
        "--percent-identity",
        type=float,
        default=97,
        help="Minimum alignment identity included in the depth summary",
    )
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ used for depth calculation")
    parser.add_argument("--weight-mapq", type=float, default=0.0, help="MAPQ weighting used for depth calculation")
    parser.add_argument("--include-edge-bases", action="store_true", help="Include contig-edge bases in depth statistics")
    parser.add_argument("--max-edge-bases", type=int, default=75, help="Edge bases excluded when edge inclusion is disabled")
    parser.add_argument("--min-contig-length", type=int, default=0, help="Minimum contig length reported")
    parser.add_argument("--min-contig-depth", type=float, default=0, help="Minimum contig depth reported")
    parser.add_argument("--bbmap-extra-args", default="", help="Additional native BBMap arguments")
    parser.add_argument("--tmp-dir", type=Path, help="Temporary directory root")
    parser.add_argument("--keep-sam", action="store_true", help="Keep the intermediate SAM file")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary mapping files")
    parser.add_argument("--print-defaults", action="store_true", help="Print machine-readable defaults and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.print_defaults:
        print("threads=16\nmemory=25g\npercent_identity=97\nmin_mapq=0\nweight_mapq=0.0")
        return 0
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if not 0 <= args.percent_identity <= 100:
        raise ValueError("--percent-identity must be between 0 and 100")
    if min(
        args.min_mapq,
        args.weight_mapq,
        args.max_edge_bases,
        args.min_contig_length,
        args.min_contig_depth,
    ) < 0:
        raise ValueError("coverage filtering thresholds must be non-negative")
    if args.long_reads:
        if args.reads1 or args.reads2:
            raise ValueError("--long-reads cannot be combined with -1 or -2")
        input_reads = (args.long_reads,)
    else:
        if not args.reads1 or not args.reads2:
            raise ValueError("paired short-read coverage requires both -1 and -2")
        input_reads = (args.reads1, args.reads2)
    for source in (*input_reads, args.reference):
        if not source.exists():
            raise FileNotFoundError(f"Input not found: {source}")
    bbmap = executable("bbmap.sh")
    samtools = executable("samtools")
    summarize = executable("jgi_summarize_bam_contig_depths")
    tmp_root = args.tmp_dir or Path(os.environ.get("METAHICT_TMP_ROOT") or os.environ.get("TMPDIR") or tempfile.gettempdir())
    tmp_root.mkdir(parents=True, exist_ok=True)
    clear_output(args.output)
    index_dir = args.output / "bbmap_index"
    index_dir.mkdir()
    run_tmp = Path(tempfile.mkdtemp(prefix="metahict_coverage.", dir=tmp_root))
    environment = os.environ.copy()
    environment.pop("JAVA_TOOL_OPTIONS", None)
    environment.pop("_JAVA_OPTIONS", None)
    environment["TMPDIR"] = str(run_tmp)
    sam = args.output / "SG_map.sam"
    sorted_bam = args.output / "SG_map_sorted.bam"
    try:
        read_arguments = (
            [f"in={args.long_reads}"]
            if args.long_reads
            else [f"in1={args.reads1}", f"in2={args.reads2}"]
        )
        bbmap_command = [bbmap, *read_arguments, f"ref={args.reference}",
                         f"path={index_dir}", f"out={sam}", f"threads={args.threads}", f"-Xmx{args.memory}",
                         *shlex.split(args.bbmap_extra_args)]
        run_logged(bbmap_command, args.output / "bbmap.log", env=environment)
        if not sam.is_file():
            raise RuntimeError("BBMap did not create SG_map.sam")
        view = subprocess.Popen([samtools, "view", "-bS", str(sam)], stdout=subprocess.PIPE)
        assert view.stdout is not None
        with (args.output / "samtools.log").open("w") as log:
            sort = subprocess.run([samtools, "sort", "-@", str(args.threads), "-T", str(run_tmp / "sort"),
                                   "-o", str(sorted_bam), "-"], stdin=view.stdout, stdout=log, stderr=subprocess.STDOUT)
        view.stdout.close()
        statuses = (view.wait(), sort.returncode)
        if any(status != 0 for status in statuses):
            raise RuntimeError(f"SAM-to-BAM conversion failed with statuses {statuses}")

        summarize_command = [summarize, "--outputDepth", str(args.output / "coverage.txt"),
                             "--percentIdentity", str(args.percent_identity),
                             "--minMapQual", str(args.min_mapq), "--weightMapQual", str(args.weight_mapq),
                             "--maxEdgeBases", str(args.max_edge_bases)]
        if not args.long_reads:
            summarize_command.extend(
                ["--pairedContigs", str(args.output / "pair.txt")]
            )
        if args.include_edge_bases:
            summarize_command.append("--includeEdgeBases")
        if args.min_contig_length > 0:
            summarize_command.extend(["--minContigLength", str(args.min_contig_length)])
        if args.min_contig_depth > 0:
            summarize_command.extend(["--minContigDepth", str(args.min_contig_depth)])
        summarize_command.append(str(sorted_bam))
        run_logged(summarize_command, args.output / "jgi_summarize.log")
        if not (args.output / "coverage.txt").is_file():
            raise RuntimeError("Coverage summarization did not create coverage.txt")
    finally:
        if not args.keep_temp:
            shutil.rmtree(run_tmp, ignore_errors=True)
            shutil.rmtree(index_dir, ignore_errors=True)
        if not args.keep_sam:
            sam.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
