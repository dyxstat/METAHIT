#!/usr/bin/env python3
"""Module 5 contact-matrix construction and normalization driver."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


def module_root(project: Path) -> Path:
    if (project / "5_contact" / "raw_contact.py").is_file():
        return project
    if (project / "modules" / "5_contact" / "raw_contact.py").is_file():
        return project / "modules"
    raise FileNotFoundError(f"Cannot locate the flat modules/5_contact directory from {project}")


def clear_output(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for child in directory.iterdir():
        if child.name == "module.log":
            continue
        if child.is_dir() and not child.is_symlink():
            shutil.rmtree(child)
        else:
            child.unlink()


def run(argv: list[str]) -> None:
    print("[INFO] Executing:", " ".join(argv), flush=True)
    subprocess.run(argv, check=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module contact --help",
    )
    parser.add_argument("method", choices=("raw", "normcc", "hiczin", "bin3c", "metator"), help="Contact normalization method")
    parser.add_argument("-p", "--project-path", type=Path, required=True, help="METAHICT repository or modules root")
    parser.add_argument("--bam", type=Path, required=True, help="Name-sorted Hi-C BAM")
    parser.add_argument("--fasta", type=Path, required=True, help="Assembly FASTA used for the BAM")
    parser.add_argument("--out", "--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument("--enzyme", required=True, help="Comma-separated restriction enzymes used by the Hi-C library")
    parser.add_argument("--metacc-min-signal", type=int, default=1, help="Minimum retained contact signal")
    parser.add_argument("--metacc-min-len", type=int, default=1000, help="Minimum contig length")
    parser.add_argument("--metacc-min-mapq", type=int, default=30, help="Minimum alignment MAPQ")
    parser.add_argument("--metacc-min-match", type=int, default=30, help="Minimum aligned match length")
    parser.add_argument("--spurious-contact-percent", "--thres", type=float, default=5, help="Lowest normalized-contact percentage removed")
    parser.add_argument("--coverage-file", type=Path, help="Coverage table required by hiczin and metator")
    parser.add_argument("--epsilon", type=float, default=1, help="Normalization epsilon")
    parser.add_argument("--max-iter", type=int, default=1000, help="Maximum matrix-balancing iterations")
    parser.add_argument("--tol", type=float, default=1e-6, help="Matrix-balancing convergence tolerance")
    parser.add_argument("--print-defaults", action="store_true", help="Print machine-readable defaults and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.print_defaults:
        print("metacc_min_signal=1\nmetacc_min_len=1000\nmetacc_min_mapq=30\nmetacc_min_match=30\nspurious_contact_percent=5")
        return 0
    if min(args.metacc_min_signal, args.metacc_min_len, args.metacc_min_mapq, args.metacc_min_match) < 0:
        raise ValueError("raw-contact filtering values must be non-negative")
    if not 0 < args.spurious_contact_percent < 100:
        raise ValueError("--spurious-contact-percent must be greater than 0 and less than 100")
    if args.epsilon <= 0 or args.max_iter < 1 or args.tol <= 0:
        raise ValueError("normalization epsilon, iterations, and tolerance must be positive")
    for source in (args.bam, args.fasta):
        if not source.exists():
            raise FileNotFoundError(f"Input not found: {source}")
    if args.method in {"hiczin", "metator"} and (args.coverage_file is None or not args.coverage_file.exists()):
        raise FileNotFoundError(f"--coverage-file is required for {args.method}")
    root = module_root(args.project_path)
    module = root / "5_contact"
    clear_output(args.out)
    (args.out / "tmp").mkdir()
    run([sys.executable, str(module / "raw_contact.py"), str(args.bam), args.enzyme, str(args.fasta), str(args.out),
         str(args.metacc_min_mapq), str(args.metacc_min_len), str(args.metacc_min_match), str(args.metacc_min_signal)])
    contig_file = args.out / "contig_info.csv"
    matrix = args.out / "Raw_contact_matrix.npz"
    if not contig_file.is_file() or not matrix.is_file():
        raise RuntimeError("Raw contact generation did not create its declared outputs")
    if args.method in {"hiczin", "metator"}:
        merged = args.out / "contig_info_with_coverage.csv"
        run([sys.executable, str(module / "add_coverage.py"), "--contig_info", str(contig_file),
             "--coverage", str(args.coverage_file), "--output", str(merged)])
        contig_file = merged
    command = [sys.executable, str(module / "normalization.py"), args.method,
               "--contig_file", str(contig_file), "--contact_matrix_file", str(matrix),
               "--output_path", str(args.out), "--thres", str(args.spurious_contact_percent),
               "--min_len", str(args.metacc_min_len), "--min_signal", str(args.metacc_min_signal)]
    if args.method in {"normcc", "hiczin", "metator"}:
        command.extend(["--epsilon", str(args.epsilon)])
    elif args.method == "bin3c":
        command.extend(["--epsilon", str(args.epsilon), "--max_iter", str(args.max_iter), "--tol", str(args.tol)])
    run(command)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
