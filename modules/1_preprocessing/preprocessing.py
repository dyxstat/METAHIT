#!/usr/bin/env python3
"""Module 1 read preprocessing."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


DEFAULTS = """Module 1 preprocessing defaults:
threads=8
prefix=derived from input R1 filename
dedup=false
minlen=50
trimq=10
qtrim=r
ftl=10
xmx=25g when called directly; Nextflow derives it from task memory
ftm=5
ktrim=r
k=23
mink=11
hdist=1
adapter_ref=BBTools adapters.fa
skip_pre_qc_report=false
skip_post_qc_report=false"""


DEFAULT_XMX = "25g"
DEFAULT_THREADS = 8


def output_prefix(read1: Path, requested: str | None) -> str:
    if requested:
        return requested
    name = read1.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name.rsplit("_", 1)[0] if "_" in name else name


def executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required executable is not on PATH: {name}")
    return path


def bbtools_adapter_reference(
    bbduk: str,
    requested: Path | None = None,
    environment: dict[str, str] | None = None,
) -> Path:
    """Locate adapters.fa without assuming one BBTools package layout."""
    if requested is not None:
        if requested.is_file():
            return requested
        raise FileNotFoundError(f"BBTools adapter reference not found: {requested}")

    env = os.environ if environment is None else environment
    executable_path = Path(bbduk).absolute()
    candidates: list[Path] = []
    conda_prefix = env.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(Path(conda_prefix) / "share" / "bbmap" / "resources" / "adapters.fa")
    # Bioconda links <prefix>/bin/bbduk*.sh into an implementation below
    # <prefix>/opt/. Do not resolve that link before deriving the prefix.
    candidates.append(executable_path.parent.parent / "share" / "bbmap" / "resources" / "adapters.fa")
    candidates.append(Path(sys.prefix) / "share" / "bbmap" / "resources" / "adapters.fa")
    # Also support an unlinked upstream BBTools tree, where resources/ is next
    # to bbduk.sh.
    candidates.append(executable_path.resolve().parent / "resources" / "adapters.fa")

    searched: list[Path] = []
    for candidate in candidates:
        if candidate in searched:
            continue
        searched.append(candidate)
        if candidate.is_file():
            return candidate
    locations = "\n  ".join(str(path) for path in searched)
    raise FileNotFoundError(f"BBTools adapter reference not found. Searched:\n  {locations}")


def run(argv: list[str], *, env: dict[str, str] | None = None) -> None:
    print("[INFO] Executing:", " ".join(argv), flush=True)
    subprocess.run(argv, check=True, env=env)


def remove_fastqc_archives(directory: Path) -> None:
    for archive in directory.glob("*.zip"):
        archive.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module preprocessing --help",
    )
    parser.add_argument("-p", "--project-path", required=True, help="METAHICT repository root")
    parser.add_argument("-1", "--r1", type=Path, required=True, help="Forward paired-end FASTQ")
    parser.add_argument("-2", "--r2", type=Path, required=True, help="Reverse paired-end FASTQ")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    parser.add_argument("--prefix", help="Output prefix; derived from R1 when omitted")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Threads passed to FastQC and BBTools")
    parser.add_argument("--minlen", type=int, default=50, help="Minimum read length retained")
    parser.add_argument("--trimq", type=int, default=10, help="BBDuk quality-trimming threshold")
    parser.add_argument("--qtrim", default="r", help="BBDuk quality-trimming direction")
    parser.add_argument("--ftl", type=int, default=10, help="Bases trimmed from the left end")
    parser.add_argument("--xmx", default=None, help="Java heap for every BBTools command, e.g. 25g")
    parser.add_argument("--ftm", type=int, default=5, help="BBDuk modulo trimming value")
    parser.add_argument("--ktrim", default="r", help="BBDuk adapter-trimming direction")
    parser.add_argument("--k", type=int, default=23, help="Adapter k-mer size")
    parser.add_argument("--mink", type=int, default=11, help="Minimum adapter k-mer size")
    parser.add_argument("--hdist", type=int, default=1, help="Adapter k-mer Hamming distance")
    parser.add_argument("--adapter-ref", type=Path, help="Adapter FASTA; auto-detected from BBTools when omitted")
    parser.add_argument(
        "--dedup",
        dest="dedup",
        action="store_true",
        help="Enable sequence-level duplicate removal with Clumpify",
    )
    parser.add_argument(
        "--no-dedup",
        dest="dedup",
        action="store_false",
        help="Disable sequence-level duplicate removal",
    )
    parser.set_defaults(dedup=False)
    parser.add_argument("--skip-pre-qc-report", action="store_true", help="Skip FastQC on input reads")
    parser.add_argument("--skip-post-qc-report", action="store_true", help="Skip FastQC on cleaned reads")
    parser.add_argument("--print-defaults", action="store_true", help="Print machine-readable defaults and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.print_defaults:
        print(DEFAULTS)
        return 0
    if args.threads < 1 or args.minlen < 1:
        raise ValueError("--threads and --minlen must be positive integers")
    if min(args.trimq, args.ftl, args.ftm, args.k, args.mink, args.hdist) < 0:
        raise ValueError("trimming thresholds and lengths must be non-negative")
    if args.k < 1 or args.mink < 1 or args.mink > args.k:
        raise ValueError("adapter k-mer lengths must be positive and --mink cannot exceed --k")
    for read in (args.r1, args.r2):
        if not read.is_file():
            raise FileNotFoundError(f"Input read file not found: {read}")

    bbduk = shutil.which("bbdukOld.sh") or shutil.which("bbduk.sh")
    if bbduk is None:
        raise RuntimeError("BBTools bbdukOld.sh or bbduk.sh is not on PATH")
    clumpify = executable("clumpify.sh")
    fastqc = None
    if not args.skip_pre_qc_report or not args.skip_post_qc_report:
        fastqc = executable("fastqc")

    adapter_ref = bbtools_adapter_reference(bbduk, args.adapter_ref)

    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix(args.r1, args.prefix)
    xmx = args.xmx or DEFAULT_XMX
    step1 = (out / f"step1_adptrim_{prefix}_1.fastq.gz", out / f"step1_adptrim_{prefix}_2.fastq.gz")
    step2 = (out / f"step2_qualtrim_{prefix}_1.fastq.gz", out / f"step2_qualtrim_{prefix}_2.fastq.gz")
    step3 = (out / f"step3_lefttrim_{prefix}_1.fastq.gz", out / f"step3_lefttrim_{prefix}_2.fastq.gz")
    final = (out / f"final_{prefix}_1.fastq.gz", out / f"final_{prefix}_2.fastq.gz")

    if fastqc and not args.skip_pre_qc_report:
        report = out / f"pre-QC_{prefix}_report"
        report.mkdir(exist_ok=True)
        run([fastqc, "-q", "-t", str(args.threads), "-o", str(report), str(args.r1), str(args.r2)])
        remove_fastqc_archives(report)

    common = [f"-Xmx{xmx}"]
    bbtools_environment = os.environ.copy()
    bbtools_environment.pop("JAVA_TOOL_OPTIONS", None)
    bbtools_environment.pop("_JAVA_OPTIONS", None)
    run([bbduk, *common, f"in1={args.r1}", f"in2={args.r2}", f"out1={step1[0]}", f"out2={step1[1]}",
         f"ref={adapter_ref}", f"ktrim={args.ktrim}", f"k={args.k}", f"mink={args.mink}",
         f"hdist={args.hdist}", f"minlen={args.minlen}", f"threads={args.threads}", "tpe", "tbo"], env=bbtools_environment)
    run([bbduk, *common, f"in1={step1[0]}", f"in2={step1[1]}", f"out1={step2[0]}", f"out2={step2[1]}",
         f"qtrim={args.qtrim}", f"trimq={args.trimq}", f"ftm={args.ftm}",
         f"minlen={args.minlen}", f"threads={args.threads}"], env=bbtools_environment)
    run([bbduk, *common, f"in1={step2[0]}", f"in2={step2[1]}", f"out1={step3[0]}", f"out2={step3[1]}",
         f"ftl={args.ftl}", f"threads={args.threads}"], env=bbtools_environment)

    if args.dedup:
        run([clumpify, *common, f"in1={step3[0]}", f"in2={step3[1]}", f"out1={final[0]}", f"out2={final[1]}",
             "dedupe", f"threads={args.threads}"], env=bbtools_environment)
    else:
        step3[0].replace(final[0])
        step3[1].replace(final[1])

    if fastqc and not args.skip_post_qc_report:
        report = out / f"post-QC_{prefix}_report"
        report.mkdir(exist_ok=True)
        run([fastqc, "-t", str(args.threads), "-o", str(report), str(final[0]), str(final[1])], env=bbtools_environment)
        remove_fastqc_archives(report)
    print(f"{prefix} preprocessing complete!")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
