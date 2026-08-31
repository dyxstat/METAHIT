#!/usr/bin/env python3
"""Integrate and quality-filter bin sets for the METAHICT binning process.

The selection workflow was adapted from metaWRAP's GPL-3.0 bin-refinement
module; provenance is recorded in docs/third_party.md. Binning_refiner itself is
an installed, checksum-pinned upstream dependency and is invoked through the
METAHICT adapter in ``run_binning_refiner.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, Optional, Sequence


class RefinementError(RuntimeError):
    """Expected bin-refinement failure."""


SCRIPT_DIR = Path(__file__).resolve().parent
HELPERS = SCRIPT_DIR


def run(
    command: Sequence[object],
    *,
    cwd: Optional[Path] = None,
    env: Optional[dict[str, str]] = None,
    capture: bool = False,
) -> subprocess.CompletedProcess[str]:
    argv = [str(item) for item in command]
    print(f"[RUN] {shlex.join(argv)}", flush=True)
    return subprocess.run(
        argv,
        cwd=str(cwd) if cwd else None,
        env=env,
        check=True,
        text=True,
        capture_output=capture,
    )


def replace_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def normalize_contig_headers(path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with path.open() as source, temporary.open("w") as output:
        for line in source:
            output.write(line.replace("=", "_") if line.startswith(">") else line)
    temporary.replace(path)


def prepare_bin_set(source: Path, destination: Path, minimum: int, maximum: int) -> int:
    if not source.is_dir():
        raise RefinementError(f"Input bin set is not a directory: {source}")
    replace_path(destination)
    destination.mkdir(parents=True)
    count = 0
    for input_file in sorted(path for path in source.iterdir() if path.is_file()):
        size = input_file.stat().st_size
        if not minimum < size < maximum:
            print(f"[INFO] Skipping {input_file}: size {size} is outside ({minimum}, {maximum})")
            continue
        output_file = destination / f"{input_file.stem}.fa"
        shutil.copy2(input_file, output_file)
        normalize_contig_headers(output_file)
        count += 1
    if not count:
        raise RefinementError(f"No valid bins remained after size filtering: {source}")
    print(f"[INFO] Prepared {count} bins in {destination.name}")
    return count


def checkm2_executable() -> str:
    executable = os.environ.get("METAHICT_CHECKM2_BIN") or shutil.which("checkm2")
    if not executable:
        raise RefinementError("CheckM2 was not found; run ./metahict install")
    return executable


def checkm2_environment(executable: str) -> dict[str, str]:
    """Expose CheckM2 companion programs from the executable's environment."""
    environment = os.environ.copy()
    executable_dir = str(Path(executable).expanduser().absolute().parent)
    path_entries = [
        entry
        for entry in environment.get("PATH", "").split(os.pathsep)
        if entry and entry != executable_dir
    ]
    environment["PATH"] = os.pathsep.join([executable_dir, *path_entries])
    if not shutil.which("diamond", path=environment["PATH"]):
        raise RefinementError(
            f"DIAMOND was not found in the CheckM2 environment: {executable_dir}"
        )
    return environment


def run_checkm2(input_dir: Path, output_dir: Path, label: str, args: argparse.Namespace) -> Path:
    replace_path(output_dir)
    temp_root = Path(args.tmp_dir).expanduser().resolve()
    temp_root.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f"metahict_checkm2_{label}_", dir=temp_root))
    executable = checkm2_executable()
    environment = checkm2_environment(executable)
    environment["TMPDIR"] = str(temporary)
    try:
        run(
            [
                executable, "predict", "-i", input_dir, "-o", output_dir,
                "-x", "fa", "-t", args.threads, "--tmpdir", temporary,
            ],
            env=environment,
        )
    finally:
        if args.keep_temp:
            print(f"[INFO] CheckM2 temporary directory retained: {temporary}")
        else:
            shutil.rmtree(temporary, ignore_errors=True)
    report = output_dir / "quality_report.tsv"
    if not report.is_file() or report.stat().st_size == 0:
        raise RefinementError(f"CheckM2 did not create a quality report: {report}")
    return report


def summarize_checkm2(
    report: Path,
    output: Path,
    *,
    binner: Optional[str] = None,
    source_stats: Optional[Path] = None,
    numeric_descending: bool = False,
) -> None:
    command: list[object] = [sys.executable, HELPERS / "summarize_checkm2.py", report]
    if source_stats is not None:
        command.extend(["manual", source_stats])
    elif binner is not None:
        command.append(binner)
    result = run(command, capture=True)
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        raise RefinementError(f"CheckM2 summarizer produced no output for {report}")
    header, records = lines[0], lines[1:]
    if numeric_descending:
        records.sort(key=lambda line: float(line.split("\t")[1]), reverse=True)
    else:
        records.sort()
    output.write_text("\n".join([header, *records]) + "\n")


def read_stats(path: Path) -> tuple[str, list[list[str]]]:
    with path.open() as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        raise RefinementError(f"Empty statistics file: {path}")
    return lines[0], [line.split("\t") for line in lines[1:]]


def good_bin_count(path: Path, completion: float, contamination: float) -> int:
    _, records = read_stats(path)
    return sum(
        len(record) >= 3
        and completion <= float(record[1]) <= 100
        and 0 <= float(record[2]) <= contamination
        for record in records
    )


def refine_one(workdir: Path, output_name: str, bin_sets: Sequence[str], minimum_size: int) -> tuple[str, int]:
    adapter = HELPERS / "run_binning_refiner.py"
    temporary = workdir / f"Refined_{output_name}"
    run(
        [sys.executable, adapter, "--output", temporary, "--min-size", minimum_size, *bin_sets],
        cwd=workdir,
    )
    refined = temporary / "Refined"
    if not refined.is_dir():
        raise RefinementError(f"Binning_refiner output is missing: {refined}")
    destination = workdir / f"bins{output_name}"
    replace_path(destination)
    refined.replace(destination)
    shutil.rmtree(temporary, ignore_errors=True)
    count = len(list(destination.glob("*.fa")))
    return destination.name, count


def all_bin_directories(workdir: Path) -> list[Path]:
    return sorted(
        path
        for path in workdir.iterdir()
        if path.is_dir() and path.name.startswith("bins") and not path.name.endswith(".checkm2")
    )


def rename_fasta_extensions(workdir: Path) -> None:
    for directory in all_bin_directories(workdir):
        for fasta in directory.glob("*.fasta"):
            fasta.rename(fasta.with_suffix(".fa"))
        if not list(directory.glob("*.fa")):
            print(f"[INFO] Removing empty refined bin set: {directory.name}")
            shutil.rmtree(directory)


def consolidate(workdir: Path, args: argparse.Namespace, number_of_inputs: int) -> str:
    if not args.consolidate:
        if not args.run_checkm2:
            return "binsA"
        candidates = [path for path in sorted(workdir.glob("bins*.stats")) if path.name != "binsM.stats"]
        if not candidates:
            raise RefinementError("No CheckM2 statistics are available to select a bin set")
        ranked = [(good_bin_count(path, args.completeness, args.contamination), path.stem) for path in candidates]
        return max(ranked, key=lambda item: (item[0], item[1]))[1]
    if number_of_inputs == 1:
        return "binsA"
    if not args.run_checkm2:
        raise RefinementError("Consolidation requires CheckM2; remove --skip-checkm2 or use --skip-consolidation")

    bins_m = workdir / "binsM"
    stats_m = workdir / "binsM.stats"
    replace_path(bins_m)
    replace_path(stats_m)
    shutil.copytree(workdir / "binsA", bins_m)
    shutil.copy2(workdir / "binsA.stats", stats_m)
    for stats in [path for path in sorted(workdir.glob("bins*.stats")) if path.name not in {"binsA.stats", "binsM.stats"}]:
        output = workdir / "binsM1"
        replace_path(output)
        replace_path(workdir / "binsM1.stats")
        run(
            [
                sys.executable, HELPERS / "consolidate_two_sets_of_bins.py",
                bins_m, workdir / stats.stem, stats_m, stats, output,
                args.completeness, args.contamination, args.contamination_penalty,
            ],
            cwd=workdir,
        )
        replace_path(bins_m)
        replace_path(stats_m)
        output.replace(bins_m)
        (workdir / "binsM1.stats").replace(stats_m)

    bins_o = workdir / "binsO"
    replace_path(bins_o)
    replace_path(workdir / "binsO.stats")
    if args.dereplication == "keep":
        bins_m.replace(bins_o)
        stats_m.replace(workdir / "binsO.stats")
    else:
        command: list[object] = [
            sys.executable, HELPERS / "dereplicate_contigs_in_bins.py", stats_m, bins_m, bins_o,
        ]
        if args.dereplication == "remove":
            command.append("remove")
        run(command, cwd=workdir)
    return "binsO"


def filter_final_bins(workdir: Path, args: argparse.Namespace) -> None:
    bins_o = workdir / "binsO"
    report = run_checkm2(bins_o, workdir / "binsO.checkm2", "binsO", args)
    stats_o = workdir / "binsO.stats"
    summarize_checkm2(
        report,
        stats_o,
        source_stats=workdir / "binsM.stats",
        numeric_descending=True,
    )
    header, records = read_stats(stats_o)
    kept: list[list[str]] = []
    for record in records:
        if len(record) < 3:
            continue
        passed = (
            args.completeness <= float(record[1]) <= 100
            and 0 <= float(record[2]) <= args.contamination
        )
        fasta = bins_o / f"{record[0]}.fa"
        if passed:
            kept.append(record)
        elif fasta.exists():
            fasta.unlink()
    stats_o.write_text(
        header + "\n" + "".join("\t".join(record) + "\n" for record in kept)
    )


def plot_stats(workdir: Path, args: argparse.Namespace, destination_name: str) -> None:
    stats = sorted(workdir.glob("*.stats"))
    if not stats:
        return
    run(
        [sys.executable, HELPERS / "plot_binning_results.py", args.completeness, args.contamination, *stats],
        cwd=workdir,
    )
    figures = workdir / "figures"
    figures.mkdir(exist_ok=True)
    for extension in ("png", "eps"):
        source = workdir / f"binning_results.{extension}"
        if source.exists():
            source.replace(figures / f"{destination_name}.{extension}")


def publish_final_outputs(
    output: Path,
    work_files: Path,
    best_bin_set: str,
    keep_intermediates: bool,
) -> Path:
    """Publish one final bin set and optionally retain refinement internals."""
    final_source = work_files / best_bin_set
    if not final_source.is_dir():
        fallback = work_files / "binsO"
        if fallback.is_dir():
            final_source = fallback
        else:
            raise RefinementError(f"Final bin set was not found: {best_bin_set}")

    final_bins = output / "final_bins"
    replace_path(final_bins)
    shutil.copytree(final_source, final_bins)
    final_stats_source = work_files / f"{final_source.name}.stats"
    if final_stats_source.is_file():
        shutil.copy2(final_stats_source, output / "final_bins_quality.tsv")

    if keep_intermediates:
        refinement_intermediates = output / "intermediates" / "refinement"
        refinement_intermediates.parent.mkdir(exist_ok=True)
        replace_path(refinement_intermediates)
        work_files.replace(refinement_intermediates)
    else:
        shutil.rmtree(work_files)
    return final_bins


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Integrate METAHICT bin sets")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-A", dest="bins_a", required=True)
    parser.add_argument("-B", dest="bins_b")
    parser.add_argument("-C", dest="bins_c")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=None,
        help="Deprecated compatibility option; task memory is controlled by the workflow",
    )
    parser.add_argument("-c", "--completeness", type=float, default=50)
    parser.add_argument("-x", "--contamination", type=float, default=10)
    parser.add_argument("--contamination-penalty", type=float, default=5)
    parser.add_argument("--min-input-bin-size", type=int, default=50000)
    parser.add_argument("--max-input-bin-size", type=int, default=20000000)
    parser.add_argument("--binning-refiner-min-size", type=int, default=500000)
    parser.add_argument("--tmp-dir", default=os.environ.get("METAHICT_TMP_ROOT", os.environ.get("TMPDIR", "/tmp")))
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep refinement and CheckM2 intermediate files for debugging",
    )
    parser.add_argument("--skip-checkm2", dest="run_checkm2", action="store_false")
    parser.add_argument("--skip-refinement", dest="refine", action="store_false")
    parser.add_argument("--skip-consolidation", dest="consolidate", action="store_false")
    dereplication = parser.add_mutually_exclusive_group()
    dereplication.add_argument("--keep-ambiguous", dest="dereplication", action="store_const", const="keep")
    dereplication.add_argument("--remove-ambiguous", dest="dereplication", action="store_const", const="remove")
    parser.set_defaults(run_checkm2=True, refine=True, consolidate=True, dereplication="best")
    return parser


def command_main(args: argparse.Namespace) -> Path:
    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    for pattern in ("binsA", "binsB", "binsC", "binsAB", "binsAC", "binsBC", "binsABC"):
        replace_path(output / pattern)
    input_paths = [Path(args.bins_a).expanduser().resolve()]
    input_paths.extend(
        Path(value).expanduser().resolve()
        for value in (args.bins_b, args.bins_c)
        if value
    )
    labels = ("binsA", "binsB", "binsC")
    for source, label in zip(input_paths, labels):
        prepare_bin_set(source, output / label, args.min_input_bin_size, args.max_input_bin_size)
    number_of_inputs = len(input_paths)

    if args.refine and number_of_inputs >= 2:
        combinations = [("AB", ("binsA", "binsB"))]
        if number_of_inputs == 3:
            combinations.extend(
                [("BC", ("binsC", "binsB")), ("AC", ("binsA", "binsC")), ("ABC", ("binsA", "binsB", "binsC"))]
            )
        with ThreadPoolExecutor(max_workers=len(combinations)) as pool:
            futures = [
                pool.submit(refine_one, output, label, bin_sets, args.binning_refiner_min_size)
                for label, bin_sets in combinations
            ]
            for future in futures:
                name, count = future.result()
                print(f"[INFO] {name}: {count} refined bins")
    rename_fasta_extensions(output)

    if args.run_checkm2:
        for directory in all_bin_directories(output):
            report = run_checkm2(directory, output / f"{directory.name}.checkm2", directory.name, args)
            stats = output / f"{directory.name}.stats"
            summarize_checkm2(report, stats, binner=directory.name)
            print(
                f"[INFO] {directory.name}: {good_bin_count(stats, args.completeness, args.contamination)} bins pass quality thresholds"
            )

    best_bin_set = consolidate(output, args, number_of_inputs)
    if args.run_checkm2 and args.consolidate and number_of_inputs > 1 and args.dereplication != "keep":
        filter_final_bins(output, args)
    if args.run_checkm2:
        plot_stats(output, args, "binning_results")

    work_files = output / "work_files"
    replace_path(work_files)
    work_files.mkdir()
    for path in sorted(output.glob("bins*")):
        if path == work_files or path.name.startswith("binning_results"):
            continue
        path.replace(work_files / path.name)
    final_bins = publish_final_outputs(
        output,
        work_files,
        best_bin_set,
        args.keep_temp,
    )
    print(f"[PASS] Final METAHICT bins: {final_bins}")
    return final_bins


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        command_main(build_parser().parse_args(argv))
    except (RefinementError, OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
