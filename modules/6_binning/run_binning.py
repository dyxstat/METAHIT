#!/usr/bin/env python3
"""Nextflow-facing entry point for METAHICT binning and bin integration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Optional, Sequence


class BinningError(RuntimeError):
    pass


DEFAULT_THREADS = 16
METAHICT_VERSION = "1.1.0"


def run(command: Sequence[object]) -> None:
    argv = [str(item) for item in command]
    print(f"[RUN] {shlex.join(argv)}", flush=True)
    subprocess.run(argv, check=True)


def write_run_parameters(destination: Path, args: argparse.Namespace) -> None:
    """Record the effective binning arguments in a human-readable manifest."""
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
        "module: binning",
        "parameters:",
    ]
    lines.extend(
        f"  {key}: {yaml_scalar(value)}"
        for key, value in vars(args).items()
    )
    destination.write_text("\n".join(lines) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run MetaCC, bin3C, and ImputeCC, then consolidate their bins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module binning --help",
    )
    parser.add_argument("--fasta", required=True, help="Assembly FASTA")
    parser.add_argument("--bam", required=True, help="Name-sorted Hi-C BAM aligned to the assembly")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--modules-path", required=True, help="METAHICT modules directory")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Threads propagated to binners, CheckM2, and refinement")
    parser.add_argument("--checkm_db", help="CheckM database directory")
    parser.add_argument("--enzyme", default="Sau3AI,MluCI", help="Comma-separated restriction enzymes")
    parser.add_argument("--metacc-min-len", type=int, help="Minimum contig length for MetaCC")
    parser.add_argument("--metacc-min-signal", type=int, help="Minimum contact signal for MetaCC")
    parser.add_argument("--metacc-min-mapq", type=int, help="Minimum MAPQ for MetaCC")
    parser.add_argument("--metacc-min-match", type=int, help="Minimum aligned match length for MetaCC")
    parser.add_argument("--metacc-min-binsize", type=int, help="Minimum output-bin size for MetaCC")
    parser.add_argument("--normcc-thres", type=float, help="NormCC denoising threshold")
    parser.add_argument("--bin3c-min-len", type=int, help="Minimum contig length for bin3C")
    parser.add_argument("--bin3c-min-signal", type=int, help="Minimum contact signal for bin3C")
    parser.add_argument("--bin3c-min-mapq", type=int, help="Minimum MAPQ for bin3C")
    parser.add_argument("--bin3c-min-match", type=int, help="Minimum aligned match length for bin3C")
    parser.add_argument("--bin3c-min-extent", type=int, help="Minimum bin3C cluster extent")
    parser.add_argument("--min-completeness", type=float, default=50, help="Minimum final-bin completeness")
    parser.add_argument("--max-contamination", type=float, default=10, help="Maximum final-bin contamination")
    parser.add_argument("--contamination-penalty", type=float, default=5, help="Contamination penalty in representative scoring")
    parser.add_argument("--min-input-bin-size", type=int, default=50000, help="Minimum candidate-bin FASTA size")
    parser.add_argument("--max-input-bin-size", type=int, default=20000000, help="Maximum candidate-bin FASTA size")
    parser.add_argument("--binning-refiner-min-size", type=int, default=500000, help="Minimum size used by Binning_refiner")
    parser.add_argument("--num-gene", type=int, help="MetaCC marker-gene count; auto-detected when omitted")
    parser.add_argument("--heatmap-max-image", type=int, default=5000, help="Maximum binning heatmap dimension")
    parser.add_argument("--imputecc-gene-coverage", type=float, default=0.9, help="ImputeCC marker-gene coverage threshold")
    parser.add_argument("--imputecc-rwr-restart-probability", type=float, default=0.5, help="ImputeCC random-walk restart probability")
    parser.add_argument("--imputecc-rwr-threshold", type=float, default=80, help="ImputeCC random-walk threshold")
    parser.add_argument("--imputecc-max-markers", type=int, default=8000, help="ImputeCC maximum marker count")
    parser.add_argument("--imputecc-intra-bin-threshold", type=float, default=50, help="ImputeCC intra-bin threshold")
    parser.add_argument("--imputecc-inter-bin-threshold", type=float, default=0, help="ImputeCC inter-bin threshold")
    parser.add_argument("--imputecc-min-bin-size", type=int, default=100000, help="ImputeCC minimum bin size")
    parser.add_argument("--imputecc-contamination-weight", type=float, default=2, help="ImputeCC contamination weight")
    parser.add_argument("--imputecc-min-completeness", type=float, default=50, help="ImputeCC minimum completeness")
    parser.add_argument("--imputecc-max-contamination", type=float, default=10, help="ImputeCC maximum contamination")
    parser.add_argument("--imputecc-report-quality-threshold", type=float, default=10, help="ImputeCC report-quality threshold")
    parser.add_argument("--tmp-dir", help="Temporary directory root")
    parser.add_argument("--seed", type=int, help="Random seed for reproducible clustering")
    boolean_help = {
        "no-fasta": "Do not write bin3C cluster FASTAs",
        "no-report": "Do not write the bin3C cluster report",
        "no-spades": "Treat input contigs as non-SPAdes identifiers",
        "only-large": "Write only bin3C clusters above the extent threshold",
        "keep-temp": "Keep refinement and CheckM2 intermediate files for debugging",
        "skip-checkm2": "Skip CheckM2 quality evaluation",
        "skip-refinement": "Skip hybrid bin refinement",
        "skip-consolidation": "Skip final non-redundant consolidation",
        "keep-ambiguous": "Keep ambiguous contigs in every candidate bin",
        "remove-ambiguous": "Remove ambiguous contigs from all candidate bins",
    }
    for option, help_text in boolean_help.items():
        parser.add_argument(f"--{option}", action="store_true", help=help_text)
    return parser


def command_main(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise BinningError("--threads must be at least 1")
    if not 0 <= args.min_completeness <= 100 or not 0 <= args.max_contamination <= 100:
        raise BinningError("completeness and contamination thresholds must be between 0 and 100")
    if args.min_input_bin_size < 0 or args.max_input_bin_size < args.min_input_bin_size:
        raise BinningError("input-bin sizes must be non-negative and maximum must be at least minimum")
    if args.binning_refiner_min_size < 0 or args.heatmap_max_image < 1:
        raise BinningError("refiner minimum size must be non-negative and heatmap size must be positive")
    if args.normcc_thres is not None and not 0 <= args.normcc_thres <= 1:
        raise BinningError("--normcc-thres must be between 0 and 1")
    optional_nonnegative = (
        args.metacc_min_len,
        args.metacc_min_signal,
        args.metacc_min_mapq,
        args.metacc_min_match,
        args.metacc_min_binsize,
        args.bin3c_min_len,
        args.bin3c_min_signal,
        args.bin3c_min_mapq,
        args.bin3c_min_match,
        args.bin3c_min_extent,
    )
    if any(value is not None and value < 0 for value in optional_nonnegative):
        raise BinningError("MetaCC and bin3C thresholds must be non-negative")
    if not 0 < args.imputecc_gene_coverage <= 1:
        raise BinningError("--imputecc-gene-coverage must be greater than 0 and at most 1")
    if not 0 <= args.imputecc_rwr_restart_probability <= 1:
        raise BinningError("--imputecc-rwr-restart-probability must be between 0 and 1")
    nonnegative = {
        "--imputecc-rwr-threshold": args.imputecc_rwr_threshold,
        "--imputecc-intra-bin-threshold": args.imputecc_intra_bin_threshold,
        "--imputecc-inter-bin-threshold": args.imputecc_inter_bin_threshold,
        "--imputecc-min-bin-size": args.imputecc_min_bin_size,
        "--imputecc-contamination-weight": args.imputecc_contamination_weight,
        "--imputecc-report-quality-threshold": args.imputecc_report_quality_threshold,
    }
    invalid = [name for name, value in nonnegative.items() if value < 0]
    if invalid or args.imputecc_max_markers < 1:
        raise BinningError(
            "ImputeCC thresholds must be non-negative and --imputecc-max-markers must be positive"
        )
    if not 0 <= args.imputecc_min_completeness <= 100 or not 0 <= args.imputecc_max_contamination <= 100:
        raise BinningError("ImputeCC completeness and contamination must be between 0 and 100")
    modules = Path(args.modules_path).expanduser().resolve()
    module = modules / "6_binning"
    binning_script = module / "generate_bins.py"
    integration_script = module / "integrate_bins.py"
    for path in (binning_script, integration_script):
        if not path.is_file():
            raise BinningError(f"Required binning implementation not found: {path}")
    fasta = Path(args.fasta).expanduser().resolve()
    bam = Path(args.bam).expanduser().resolve()
    for path in (fasta, bam):
        if not path.is_file():
            raise BinningError(f"Input file not found: {path}")
    output = Path(args.outdir).expanduser().resolve()
    if output.exists():
        shutil.rmtree(output)

    binning_command: list[object] = [
        sys.executable, binning_script,
        "--FASTA", fasta, "--BAM", bam, "--OUTDIR", output,
        "--threads", args.threads, "--enzyme", args.enzyme,
    ]
    if args.checkm_db:
        binning_command.extend(["--checkm_db", Path(args.checkm_db).expanduser().resolve()])
    value_options = (
        ("metacc_min_len", "--metacc-min-len"),
        ("metacc_min_signal", "--metacc-min-signal"),
        ("metacc_min_mapq", "--metacc-min-mapq"),
        ("metacc_min_match", "--metacc-min-match"),
        ("metacc_min_binsize", "--metacc-min-binsize"),
        ("normcc_thres", "--thres"),
        ("bin3c_min_len", "--bin3c-min-len"),
        ("bin3c_min_signal", "--bin3c-min-signal"),
        ("bin3c_min_mapq", "--bin3c-min-mapq"),
        ("bin3c_min_match", "--bin3c-min-match"),
        ("bin3c_min_extent", "--bin3c-min-extent"),
        ("num_gene", "--num-gene"),
        ("imputecc_gene_coverage", "--imputecc-gene-coverage"),
        ("imputecc_rwr_restart_probability", "--imputecc-rwr-restart-probability"),
        ("imputecc_rwr_threshold", "--imputecc-rwr-threshold"),
        ("imputecc_max_markers", "--imputecc-max-markers"),
        ("imputecc_intra_bin_threshold", "--imputecc-intra-bin-threshold"),
        ("imputecc_inter_bin_threshold", "--imputecc-inter-bin-threshold"),
        ("imputecc_min_bin_size", "--imputecc-min-bin-size"),
        ("imputecc_contamination_weight", "--imputecc-contamination-weight"),
        ("imputecc_min_completeness", "--imputecc-min-completeness"),
        ("imputecc_max_contamination", "--imputecc-max-contamination"),
        ("imputecc_report_quality_threshold", "--imputecc-report-quality-threshold"),
        ("seed", "--seed"),
    )
    for attribute, option in value_options:
        value = getattr(args, attribute)
        if value is not None:
            binning_command.extend([option, value])
    for attribute, option in (
        ("no_fasta", "--no-fasta"), ("no_report", "--no-report"),
        ("no_spades", "--no-spades"), ("only_large", "--only-large"),
    ):
        if getattr(args, attribute):
            binning_command.append(option)
    run(binning_command)

    integration_command: list[object] = [
        sys.executable, integration_script,
        output / "metacc", output / "bin3c", output / "imputecc", output / "metahict",
        "--threads", args.threads,
        "--min-completeness", args.min_completeness,
        "--max-contamination", args.max_contamination,
        "--contamination-penalty", args.contamination_penalty,
        "--min-input-bin-size", args.min_input_bin_size,
        "--max-input-bin-size", args.max_input_bin_size,
        "--binning-refiner-min-size", args.binning_refiner_min_size,
    ]
    if args.tmp_dir:
        integration_command.extend(["--tmp-dir", Path(args.tmp_dir).expanduser().resolve()])
    for attribute, option in (
        ("keep_temp", "--keep-temp"), ("skip_checkm2", "--skip-checkm2"),
        ("skip_refinement", "--skip-refinement"), ("skip_consolidation", "--skip-consolidation"),
        ("keep_ambiguous", "--keep-ambiguous"), ("remove_ambiguous", "--remove-ambiguous"),
    ):
        if getattr(args, attribute):
            integration_command.append(option)
    run(integration_command)

    metahict_output = output / "metahict"
    final_bins = metahict_output / "final_bins"
    final_quality = metahict_output / "final_bins_quality.tsv"
    combined_bins = metahict_output / "combined_final_bins.fa"
    membership = metahict_output / "contig_to_bin.tsv"
    if not final_bins.is_dir() or not any(final_bins.glob("*.fa")):
        raise BinningError(f"Expected consolidated bins were not created: {final_bins}")
    for required in (combined_bins, membership):
        if not required.is_file() or required.stat().st_size == 0:
            raise BinningError(f"Expected binning result was not created: {required}")
    if not args.skip_checkm2 and (
        not final_quality.is_file() or final_quality.stat().st_size == 0
    ):
        raise BinningError(f"Expected final-bin quality table was not created: {final_quality}")
    write_run_parameters(metahict_output / "run_parameters.yaml", args)

    heatmap = module / "heatmap.py"
    figure_dir = metahict_output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    serialized_membership = metahict_output / "intermediates" / "final_bins.p.gz"
    run(
        [
            sys.executable, heatmap,
            "--contact-map", output / "metacc" / "Normalized_contact_matrix.npz",
            "--ORDER", output / "metacc" / "contact_map.p.gz",
            "--BIN", serialized_membership,
            "--OUTDIR", figure_dir,
            "--max_image", args.heatmap_max_image,
        ]
    )
    if not args.keep_temp:
        shutil.rmtree(metahict_output / "intermediates", ignore_errors=True)
    print(f"[PASS] Binning completed: {metahict_output}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        command_main(build_parser().parse_args(argv))
    except (BinningError, OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
