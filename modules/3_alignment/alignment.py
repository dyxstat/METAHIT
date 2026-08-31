#!/usr/bin/env python3
"""Align Hi-C reads to an assembly and calculate alignment quality metrics."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
DEFAULT_THREADS = 16


def executable(name: str) -> str:
    result = shutil.which(name)
    if result is None:
        raise RuntimeError(f"Required executable is not on PATH: {name}")
    return result


def match_length(cigar: str) -> int:
    if cigar == "*":
        return 0
    return sum(int(length) for length, operation in CIGAR_RE.findall(cigar) if operation in {"M", "=", "X"})


def informative_pair(first: list[str], second: list[str], min_distance: int, min_match: int, *, check_distance: bool) -> bool:
    flag_a, flag_b = int(first[1]), int(second[1])
    if flag_a & (0x4 | 0x400) or flag_b & (0x4 | 0x400):
        return False
    if int(first[4]) == 0 or int(second[4]) == 0:
        return False
    if match_length(first[5]) < min_match or match_length(second[5]) < min_match:
        return False
    if first[2] != second[2]:
        return True
    return check_distance and abs(int(first[3]) - int(second[3])) >= min_distance


def write_metrics(samtools: str, bam: Path, output: Path, mapq: int, min_distance: int, min_match: int) -> None:
    process = subprocess.Popen(
        [samtools, "view", "-f", "1", "-F", "0x904", "-q", str(mapq), str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    total = informative = chimeric = 0
    first: list[str] | None = None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        if first is None:
            first = fields
        elif fields[0] == first[0]:
            total += 1
            informative += informative_pair(first, fields, min_distance, min_match, check_distance=True)
            chimeric += informative_pair(first, fields, min_distance, min_match, check_distance=False)
            first = None
        else:
            first = fields
    if process.wait() != 0:
        raise RuntimeError("samtools failed while calculating alignment metrics")
    ratio_3d = float("nan") if total == 0 else chimeric / total
    informative_ratio = 0.0 if total == 0 else informative / total
    (output / "3d_ratio.txt").write_text(f"{ratio_3d:.4f}")
    (output / "informative_pairs_ratio.txt").write_text(f"{informative_ratio:.4f}")
    print(f"[INFO] Total pairs: {total}; informative: {informative}; chimeric: {chimeric}")


def filter_alignment_stream(input_handle, output_handle, min_match: int) -> None:
    for line in input_handle:
        if line.startswith("@"):
            output_handle.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 6 and (min_match <= 0 or match_length(fields[5]) >= min_match):
            output_handle.write(line)


def sam_to_sorted_bam(
    samtools: str,
    sam_source,
    output: Path,
    filters: list[str],
    mapq: int,
    min_match: int,
    threads: int,
    sort_memory: str,
    temp_prefix: Path,
) -> None:
    view = subprocess.Popen(
        [samtools, "view", *filters, "-h", "-q", str(mapq), sam_source],
        stdout=subprocess.PIPE,
        text=True,
    )
    binary = subprocess.Popen([samtools, "view", "-bS", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    assert view.stdout is not None and binary.stdin is not None and binary.stdout is not None
    sort = subprocess.Popen(
        [samtools, "sort", "-n", "-@", str(threads), "-m", sort_memory, "-T", str(temp_prefix), "-o", str(output), "-"],
        stdin=binary.stdout,
    )
    binary.stdout.close()
    try:
        filter_alignment_stream(view.stdout, binary.stdin, min_match)
    finally:
        binary.stdin.close()
    statuses = (view.wait(), binary.wait(), sort.wait())
    if any(status != 0 for status in statuses):
        raise RuntimeError(f"samtools alignment pipeline failed with statuses {statuses}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module alignment --help",
    )
    parser.add_argument("-p", "--project-path", required=True, help="METAHICT repository root")
    parser.add_argument("-r", "--reference", type=Path, required=True, help="Assembly FASTA")
    parser.add_argument("-1", "--reads1", type=Path, required=True, help="Cleaned Hi-C forward FASTQ")
    parser.add_argument("-2", "--reads2", type=Path, required=True, help="Cleaned Hi-C reverse FASTQ")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Threads passed to BWA-MEM and SAMtools")
    parser.add_argument("--bwa-options", default="-5SP", help="Options passed to BWA-MEM")
    parser.add_argument("--samtools-filter", default="-F 0x900", help="Options passed to samtools view")
    parser.add_argument("--mapq", type=int, default=30, help="Minimum mapping quality retained")
    parser.add_argument("--min-intra-dist", type=int, default=10000, help="Minimum intra-contig distance used by metrics")
    parser.add_argument("--min-match-len", type=int, default=30, help="Minimum aligned match length retained")
    parser.add_argument("--sort-memory", default="1G", help="Memory per SAMtools sort thread")
    parser.add_argument("--tmp-dir", type=Path, help="Temporary directory root")
    parser.add_argument("--keep-sam", action="store_true", help="Keep the intermediate SAM file")
    parser.add_argument("--skip-metrics", action="store_true", help="Skip 3D and informative-pair metrics")
    parser.add_argument("--print-defaults", action="store_true", help="Print machine-readable defaults and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.print_defaults:
        print("threads=16\nbwa_options=-5SP\nsamtools_filter=-F 0x900\nmapq=30\nmin_intra_dist=10000\nmin_match_len=30\nsort_memory=1G")
        return 0
    if min(args.threads, args.mapq, args.min_intra_dist, args.min_match_len) < 0 or args.threads < 1:
        raise ValueError("thread and alignment filter values must be non-negative; threads must be positive")
    for source in (args.reference, args.reads1, args.reads2):
        if not source.is_file():
            raise FileNotFoundError(f"Input file not found: {source}")
    bwa, samtools = executable("bwa"), executable("samtools")
    tmp_root = args.tmp_dir or Path(os.environ.get("METAHICT_TMP_ROOT") or os.environ.get("TMPDIR") or tempfile.gettempdir())
    if not tmp_root.is_dir() or not os.access(tmp_root, os.W_OK):
        raise RuntimeError(f"Temporary directory is not writable: {tmp_root}")
    args.output.mkdir(parents=True, exist_ok=True)
    reference_dir = args.output / "reference"
    reference_dir.mkdir(exist_ok=True)
    work_reference = reference_dir / args.reference.name
    shutil.copy2(args.reference, work_reference)
    run_tmp = Path(tempfile.mkdtemp(prefix="metahict_alignment.", dir=tmp_root))
    try:
        subprocess.run([bwa, "index", str(work_reference)], check=True)
        bwa_command = [bwa, "mem", *shlex.split(args.bwa_options), "-t", str(args.threads), str(work_reference), str(args.reads1), str(args.reads2)]
        map_sam = args.output / "map.sam"
        if args.keep_sam:
            with map_sam.open("w") as handle:
                subprocess.run(bwa_command, check=True, stdout=handle)
            sam_source = str(map_sam)
            sam_to_sorted_bam(samtools, sam_source, args.output / "sorted_map.bam", shlex.split(args.samtools_filter),
                              args.mapq, args.min_match_len, args.threads, args.sort_memory, run_tmp / "sort")
        else:
            bwa_process = subprocess.Popen(bwa_command, stdout=subprocess.PIPE, text=True)
            assert bwa_process.stdout is not None
            # Feed BWA directly into the first samtools process while keeping argument handling explicit.
            view = subprocess.Popen([samtools, "view", *shlex.split(args.samtools_filter), "-h", "-q", str(args.mapq), "-"],
                                    stdin=bwa_process.stdout, stdout=subprocess.PIPE, text=True)
            bwa_process.stdout.close()
            binary = subprocess.Popen([samtools, "view", "-bS", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            assert view.stdout is not None and binary.stdin is not None and binary.stdout is not None
            sort = subprocess.Popen([samtools, "sort", "-n", "-@", str(args.threads), "-m", args.sort_memory,
                                     "-T", str(run_tmp / "sort"), "-o", str(args.output / "sorted_map.bam"), "-"], stdin=binary.stdout)
            binary.stdout.close()
            try:
                filter_alignment_stream(view.stdout, binary.stdin, args.min_match_len)
            finally:
                binary.stdin.close()
            statuses = (bwa_process.wait(), view.wait(), binary.wait(), sort.wait())
            if any(status != 0 for status in statuses):
                raise RuntimeError(f"alignment pipeline failed with statuses {statuses}")
        bam = args.output / "sorted_map.bam"
        if not bam.is_file() or bam.stat().st_size == 0:
            raise RuntimeError("Alignment did not create a non-empty sorted_map.bam")
        if not args.skip_metrics:
            write_metrics(samtools, bam, args.output, args.mapq, args.min_intra_dist, args.min_match_len)
    finally:
        shutil.rmtree(run_tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
