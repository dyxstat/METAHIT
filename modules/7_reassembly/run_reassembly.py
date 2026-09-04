#!/usr/bin/env python3
"""Nextflow-facing entry point for METAHICT reassembly."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Optional, Sequence


class ReassemblyLaunchError(RuntimeError):
    pass


DEFAULT_THREADS = 16
DEFAULT_MEMORY_GB = 51
DEFAULT_EM_INITIAL_N_FRACTION = 0.8
DEFAULT_EM_CONVERGENCE_TOLERANCE = 0.01
DEFAULT_EM_MAX_ITERATIONS = 100
METAHICT_VERSION = "1.2.0"


def run(command: Sequence[object]) -> None:
    argv = [str(item) for item in command]
    print(f"[RUN] {shlex.join(argv)}", flush=True)
    subprocess.run(argv, check=True)


def replace_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def write_run_parameters(destination: Path, args: argparse.Namespace) -> None:
    def yaml_scalar(value: object) -> str:
        if value is None:
            return "null"
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (int, float)):
            return str(value)
        return json.dumps(str(value), ensure_ascii=False)

    lines = [
        f"metahict_version: {yaml_scalar(METAHICT_VERSION)}",
        "module: reassembly",
        "parameters:",
    ]
    lines.extend(
        f"  {key}: {yaml_scalar(value)}"
        for key, value in vars(args).items()
    )
    destination.write_text("\n".join(lines) + "\n")


def relocate_or_remove(
    paths: Sequence[Path],
    destination: Optional[Path],
) -> None:
    available = [
        source
        for source in paths
        if source.exists() or source.is_symlink()
    ]
    if destination is not None and available:
        destination.mkdir(parents=True, exist_ok=True)
    for source in available:
        if destination is None:
            replace_path(source)
            continue
        target = destination / source.name
        replace_path(target)
        source.replace(target)


def finalize_output_layout(output: Path, args: argparse.Namespace) -> None:
    """Keep stable results and remove or organize module-created work files."""
    final_bins = output / "reassembled_bins"
    combined = output / "combined_contigs.fa"
    summary = output / "read_selection_summary.json"
    if not final_bins.is_dir() or not any(final_bins.glob("*.fa")):
        raise ReassemblyLaunchError(f"Reassembled bins were not created: {final_bins}")
    for required in (combined, summary):
        if not required.is_file() or required.stat().st_size == 0:
            raise ReassemblyLaunchError(f"Expected reassembly result was not created: {required}")
    if not args.skip_checkm2:
        quality = output / "reassembled_bins_quality.tsv"
        if not quality.is_file() or quality.stat().st_size == 0:
            raise ReassemblyLaunchError(f"Final-bin quality table was not created: {quality}")

    figures = output / "figures"
    figure_paths = [
        output / name
        for name in (
            "reassembled_bins.png",
            "reassembled_bins_scatter.png",
            "reassembly_results.png",
            "reassembly_results.eps",
        )
    ]
    relocate_or_remove(figure_paths, figures)

    read_selection_paths = [
        output / name
        for name in (
            "hic.name_sorted.bam",
            "all_intra_insert_sizes.tsv.gz",
            "readname_sg_in_hic.txt",
            "sg_in_hic.forward.fastq.gz",
            "sg_in_hic.reverse.fastq.gz",
            "new_sg_forward.fastq.gz",
            "new_sg_reverse.fastq.gz",
            "all_hic_readnames.txt",
            "readname_non_sg_in_hic.txt",
            "non_sg_in_hic.forward.fastq.gz",
            "non_sg_in_hic.reverse.fastq.gz",
            # Obsolete names are included so reruns cannot republish them.
            "insert_size_counts.json",
            "em_parameters.json",
            "em_selection_summary.json",
            "mixing_proportion.txt",
            "insert_size_cutoff.txt",
            "informative_fraction.txt",
            "long_range_ratio.txt",
            "3d_ratio.txt",
        )
    ]
    read_selection_paths.extend(output.glob("top*_intra_insert_sizes.tsv.gz"))
    read_selection_paths.extend(output.glob("top_*_contigs.txt"))

    reassembly_paths = [
        output / name
        for name in (
            "input_assembly",
            "contig_sort.txt",
            "original_bins",
            "binned_assembly",
            "reads_for_reassembly",
            "reassemblies",
            "reassembled_best_bins",
            "reassembled_bins.checkm2",
            "residual_assembly",
            "unmapped_shotgun_1.fastq",
            "unmapped_shotgun_2.fastq",
            "work_files",
        )
    ]
    if args.keep_temp:
        intermediates = output / "intermediates"
        relocate_or_remove(read_selection_paths, intermediates / "read_selection")
        relocate_or_remove(reassembly_paths, intermediates / "reassembly")
    else:
        relocate_or_remove(read_selection_paths, None)
        relocate_or_remove(reassembly_paths, None)
        replace_path(output / "intermediates")

    write_run_parameters(output / "run_parameters.yaml", args)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Select shotgun-like Hi-C reads and reassemble each MAG",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module reassembly --help",
    )
    parser.add_argument("-p", "--modules-path", required=True, help="METAHICT modules directory")
    parser.add_argument("--bin", required=True, help="Directory containing consolidated input-bin FASTAs")
    parser.add_argument("--assembly", required=True, help="Original assembly FASTA")
    parser.add_argument("--hic1", required=True, help="Cleaned Hi-C forward FASTQ")
    parser.add_argument("--hic2", required=True, help="Cleaned Hi-C reverse FASTQ")
    parser.add_argument("--sg1", required=True, help="Cleaned shotgun forward FASTQ")
    parser.add_argument("--sg2", required=True, help="Cleaned shotgun reverse FASTQ")
    parser.add_argument("--bam", required=True, help="Name-sorted Hi-C BAM aligned to the assembly")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Total threads shared by recruitment, CheckM2, and SPAdes")
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=DEFAULT_MEMORY_GB,
        help="Total SPAdes memory budget in GB, divided among parallel bin assemblies",
    )
    parser.add_argument("--cutoff-quantile", type=float, default=0.95, help="EM short-insert cutoff quantile")
    parser.add_argument("--top-k", type=int, default=100, help="Longest contigs used for EM fitting")
    parser.add_argument(
        "--em-initial-n-fraction",
        type=float,
        default=DEFAULT_EM_INITIAL_N_FRACTION,
        help="Fraction of insert sizes used to initialize the lower N component",
    )
    parser.add_argument(
        "--em-convergence-tolerance",
        type=float,
        default=DEFAULT_EM_CONVERGENCE_TOLERANCE,
        help="Absolute log-likelihood change required for EM convergence",
    )
    parser.add_argument(
        "--em-max-iterations",
        type=int,
        default=DEFAULT_EM_MAX_ITERATIONS,
        help="Maximum number of EM fitting iterations",
    )
    parser.add_argument("--min-mapq", type=int, default=30, help="Minimum MAPQ for insert-size extraction")
    parser.add_argument("--min-match-len", type=int, default=30, help="Minimum aligned match length")
    parser.add_argument("--exclude-duplicates", action="store_true", help="Exclude duplicate-marked alignments")
    parser.add_argument("--write-nonselected-hic", action="store_true", help="Write Hi-C reads rejected by the EM selection")
    parser.add_argument("--min-contig-len", type=int, default=500, help="Minimum reassembled contig length")
    parser.add_argument("--strict-cut-off", type=int, default=2, help="Strict mismatch cutoff for read recruitment")
    parser.add_argument("--permissive-cut-off", type=int, default=5, help="Permissive mismatch cutoff for read recruitment")
    parser.add_argument("--contamination-penalty", type=float, default=5, help="Contamination penalty in representative scoring")
    parser.add_argument("--min-completeness", type=float, default=50, help="Minimum completeness for accepting a reassembled bin")
    parser.add_argument("--max-contamination", type=float, default=10, help="Maximum contamination for accepting a reassembled bin")
    parser.add_argument("--skip-checkm2", action="store_true", help="Skip CheckM2 quality evaluation")
    parser.add_argument("--checkm2_db", help="CheckM2 DIAMOND database file")
    parser.add_argument("--tmp-dir", help="Temporary directory root")
    parser.add_argument("--spades-mode", choices=("careful", "none"), default="careful", help="SPAdes mismatch-correction mode")
    parser.add_argument("--spades-phred-offset", help="Explicit SPAdes PHRED offset")
    parser.add_argument("--spades-extra-args", default="", help="Additional native SPAdes arguments")
    parser.add_argument("--skip-residual-assembly", action="store_true", help="Do not assemble residual unassigned reads")
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep read-selection, reassembly, SPAdes, and CheckM2 intermediate files",
    )
    return parser


def command_main(args: argparse.Namespace) -> None:
    if args.threads < 1 or args.memory < 1:
        raise ReassemblyLaunchError("--threads and --memory must be at least 1")
    if not 0 < args.cutoff_quantile < 1:
        raise ReassemblyLaunchError("--cutoff-quantile must be greater than 0 and less than 1")
    if not 0 < args.em_initial_n_fraction < 1:
        raise ReassemblyLaunchError(
            "--em-initial-n-fraction must be greater than 0 and less than 1"
        )
    if args.em_convergence_tolerance <= 0:
        raise ReassemblyLaunchError("--em-convergence-tolerance must be greater than 0")
    if args.em_max_iterations < 1:
        raise ReassemblyLaunchError("--em-max-iterations must be at least 1")
    if not 0 <= args.min_completeness <= 100 or not 0 <= args.max_contamination <= 100:
        raise ReassemblyLaunchError("completeness and contamination thresholds must be between 0 and 100")
    modules_path = Path(args.modules_path).expanduser().resolve()
    runner = modules_path / "7_reassembly" / "select_reassembly_reads.py"
    if not runner.is_file():
        raise ReassemblyLaunchError(f"Reassembly implementation not found: {runner}")
    output = Path(args.outdir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    assembly = Path(args.assembly).expanduser().resolve()
    if not assembly.is_file():
        raise ReassemblyLaunchError(f"Assembly not found: {assembly}")
    workdir = output / "input_assembly"
    workdir.mkdir(exist_ok=True)
    working_assembly = workdir / assembly.name
    shutil.copy2(assembly, working_assembly)
    samtools = shutil.which("samtools")
    if not samtools:
        raise ReassemblyLaunchError("samtools was not found in PATH")
    index = Path(f"{working_assembly}.fai")
    if not index.is_file():
        run([samtools, "faidx", working_assembly])
    records = []
    with index.open() as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 2:
                records.append((fields[0], int(fields[1])))
    records.sort(key=lambda item: (-item[1], item[0]))
    contig_sort = output / "contig_sort.txt"
    contig_sort.write_text("".join(f"{name}\n" for name, _ in records))
    command: list[object] = [
        sys.executable, runner,
        "--bin", Path(args.bin).expanduser().resolve(),
        "--hic1", Path(args.hic1).expanduser().resolve(),
        "--hic2", Path(args.hic2).expanduser().resolve(),
        "--sg1", Path(args.sg1).expanduser().resolve(),
        "--sg2", Path(args.sg2).expanduser().resolve(),
        "--bam", Path(args.bam).expanduser().resolve(),
        "--outdir", output,
        "-p", modules_path,
        "-c", contig_sort,
        "-t", args.threads,
        "-m", args.memory,
        "--top_k", args.top_k,
        "--min-mapq", args.min_mapq,
        "--min-match-len", args.min_match_len,
        "--cutoff-quantile", args.cutoff_quantile,
        "--em-initial-n-fraction", args.em_initial_n_fraction,
        "--em-convergence-tolerance", args.em_convergence_tolerance,
        "--em-max-iterations", args.em_max_iterations,
        "--min-contig-len", args.min_contig_len,
        "--strict-cut-off", args.strict_cut_off,
        "--permissive-cut-off", args.permissive_cut_off,
        "--contamination-penalty", args.contamination_penalty,
        "--min-completeness", args.min_completeness,
        "--max-contamination", args.max_contamination,
        "--spades-mode", args.spades_mode,
        "--bam-name-sorted",
    ]
    for enabled, option in (
        (args.exclude_duplicates, "--exclude-duplicates"),
        (args.write_nonselected_hic, "--write-nonselected-hic"),
        (args.skip_checkm2, "--skip-checkm2"),
        (args.skip_residual_assembly, "--skip-residual-assembly"),
        (args.keep_temp, "--keep-temp"),
    ):
        if enabled:
            command.append(option)
    for value, option in (
        (args.checkm2_db, "--checkm2_db"),
        (args.tmp_dir, "--tmp-dir"),
        (args.spades_phred_offset, "--spades-phred-offset"),
        (args.spades_extra_args, "--spades-extra-args"),
    ):
        if value:
            command.extend([option, value])
    run(command)
    finalize_output_layout(output, args)
    print(f"[PASS] Reassembly completed: {output}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        command_main(build_parser().parse_args(argv))
    except (ReassemblyLaunchError, OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
