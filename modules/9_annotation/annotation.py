#!/usr/bin/env python3
"""Classify metagenome-assembled genomes (MAGs) with GTDB-Tk."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


DEFAULTS = """Module 9 annotation defaults:
threads=8
pplacer_cpus=same as threads
gtdbtk_db=<project_path>/databases/gtdbtk_db/release220
extension=fa
prefix=gtdbtk
skip_ani_screen=true
mash_db=
min_perc_aa=10
min_af=0.5
full_tree=false
scratch_dir=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
force=false
keep_intermediates=false
debug=false
write_single_copy_genes=false
gtdbtk_extra_args="""


def find_gtdbtk(project: Path) -> str:
    installed = shutil.which("gtdbtk")
    if installed:
        return installed
    bundled_environment = project / "conda_envs" / "gtdbtk-2.4.0" / "bin" / "gtdbtk"
    if bundled_environment.is_file() and os.access(bundled_environment, os.X_OK):
        return str(bundled_environment)
    raise RuntimeError("gtdbtk is not on PATH; install the pinned METAHICT environments")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module annotation --help",
    )
    parser.add_argument("-p", "--project-path", type=Path, required=True, help="METAHICT repository root")
    parser.add_argument(
        "--mag-dir",
        type=Path,
        dest="mag_dir",
        required=True,
        help="Directory containing MAG FASTA files",
    )
    parser.add_argument("--outdir", type=Path, required=True, help="GTDB-Tk output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Threads passed to GTDB-Tk as --cpus")
    parser.add_argument("--pplacer-cpus", type=int, help="Threads passed specifically to pplacer; defaults to --threads")
    parser.add_argument("--gtdbtk_db", "--gtdbtk-db", type=Path, help="GTDB-Tk reference-data directory")
    parser.add_argument("--extension", default="fa", help="Genome FASTA filename extension")
    parser.add_argument("--prefix", default="gtdbtk", help="GTDB-Tk output prefix")
    parser.add_argument("--skip-ani-screen", action="store_true", default=True, help="Skip ANI/Mash prescreening")
    parser.add_argument("--no-skip-ani-screen", dest="skip_ani_screen", action="store_false", help="Enable ANI/Mash prescreening")
    parser.add_argument("--mash-db", type=Path, help="Mash database required when ANI screening is enabled")
    parser.add_argument("--min-perc-aa", type=float, default=10, help="Minimum amino-acid percentage in the MSA")
    parser.add_argument("--min-af", type=float, default=0.5, help="Minimum alignment fraction for species assignment")
    parser.add_argument("--full-tree", action="store_true", help="Use the full bacterial reference tree")
    parser.add_argument("--scratch-dir", type=Path, help="Disk-backed pplacer scratch directory")
    parser.add_argument("--tmp-dir", type=Path, help="Temporary directory root")
    parser.add_argument("--force", action="store_true", help="Continue when one genome fails")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep GTDB-Tk intermediate files")
    parser.add_argument("--debug", action="store_true", help="Enable GTDB-Tk debug mode")
    parser.add_argument("--write-single-copy-genes", action="store_true", help="Write unaligned single-copy marker genes")
    parser.add_argument("--gtdbtk-extra-args", default="", help="Additional native gtdbtk classify_wf arguments")
    parser.add_argument("--print-defaults", action="store_true", help="Print machine-readable defaults and exit")
    return parser


def build_command(args: argparse.Namespace, executable: str, database: Path, tmp_dir: Path) -> list[str]:
    command = [executable, "classify_wf", "--genome_dir", str(args.mag_dir), "--out_dir", str(args.outdir),
               "--extension", args.extension, "--prefix", args.prefix, "--cpus", str(args.threads),
               "--pplacer_cpus", str(args.pplacer_cpus or args.threads), "--min_perc_aa", str(args.min_perc_aa),
               "--min_af", str(args.min_af), "--tmpdir", str(tmp_dir)]
    if args.skip_ani_screen:
        command.append("--skip_ani_screen")
    else:
        if args.mash_db is None:
            raise ValueError("--mash-db is required with --no-skip-ani-screen")
        command.extend(["--mash_db", str(args.mash_db)])
    if args.full_tree:
        command.append("--full_tree")
    if args.scratch_dir:
        args.scratch_dir.mkdir(parents=True, exist_ok=True)
        command.extend(["--scratch_dir", str(args.scratch_dir)])
    for enabled, option in (
        (args.force, "--force"),
        (args.keep_intermediates, "--keep_intermediates"),
        (args.debug, "--debug"),
        (args.write_single_copy_genes, "--write_single_copy_genes"),
    ):
        if enabled:
            command.append(option)
    command.extend(shlex.split(args.gtdbtk_extra_args))
    return command


def main(argv: list[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    if "--print-defaults" in raw_args:
        print(DEFAULTS)
        return 0
    args = build_parser().parse_args(raw_args)
    if args.mash_db is not None:
        args.skip_ani_screen = False
    if args.threads < 1 or (args.pplacer_cpus is not None and args.pplacer_cpus < 1):
        raise ValueError("thread counts must be positive")
    if args.pplacer_cpus is not None and args.pplacer_cpus > args.threads:
        raise ValueError("--pplacer-cpus cannot exceed --threads")
    if not 0 <= args.min_perc_aa <= 100 or not 0 <= args.min_af <= 1:
        raise ValueError("--min-perc-aa must be between 0 and 100 and --min-af between 0 and 1")
    if not args.mag_dir.is_dir():
        raise FileNotFoundError(f"MAG directory not found: {args.mag_dir}")
    extension = args.extension.lstrip(".")
    if not extension:
        raise ValueError("--extension must not be empty")
    if not any(path.is_file() for path in args.mag_dir.glob(f"*.{extension}")):
        raise FileNotFoundError(
            f"No *.{extension} genome files found in MAG directory: {args.mag_dir}"
        )
    database = args.gtdbtk_db or args.project_path / "databases" / "gtdbtk_db" / "release220"
    if not database.is_dir():
        raise FileNotFoundError(f"GTDB-Tk database not found: {database}")
    tmp_dir = args.tmp_dir or Path(os.environ.get("METAHICT_TMP_ROOT") or os.environ.get("TMPDIR") or tempfile.gettempdir())
    tmp_dir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    gtdbtk = find_gtdbtk(args.project_path)
    environment = os.environ.copy()
    environment["GTDBTK_DATA_PATH"] = str(database)
    environment["TMPDIR"] = str(tmp_dir)
    command = build_command(args, gtdbtk, database, tmp_dir)
    print("[INFO] Executing:", shlex.join(command), flush=True)
    subprocess.run(command, check=True, env=environment)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
