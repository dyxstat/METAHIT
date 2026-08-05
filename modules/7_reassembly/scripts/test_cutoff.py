#!/usr/bin/env python3
import argparse
import importlib.util
import os
import shutil
import subprocess


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REASSEMBLY_SCRIPT = os.path.join(SCRIPT_DIR, "7_reassembly.py")

spec = importlib.util.spec_from_file_location("metahict_reassembly", REASSEMBLY_SCRIPT)
reassembly = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reassembly)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test SG-like Hi-C insert-size cutoffs from an existing reassembly output directory."
    )
    parser.add_argument("--outdir", required=True, help="Existing reassembly output directory containing reassembly cache files")
    parser.add_argument("--bam", required=True, help="Hi-C BAM mapped to the assembly")
    parser.add_argument("--hic1", help="Hi-C R1 FASTQ, required only with --write-fastq or missing cached counts")
    parser.add_argument("--hic2", help="Hi-C R2 FASTQ, required only with --write-fastq")
    parser.add_argument("--sg1", help="Shotgun R1 FASTQ, required only with --write-fastq or missing cached counts")
    parser.add_argument("--sg2", help="Shotgun R2 FASTQ, required only with --write-fastq")
    parser.add_argument(
        "--thresholds",
        default="50,55,60,65,70,75,80,85,90",
        help="Comma-separated percentile thresholds to test (default: 50,55,...,90)",
    )
    parser.add_argument(
        "--write-fastq",
        action="store_true",
        help="Also write threshold-specific sg_in_hic/new_sg/non_sg FASTQ files. This can be large.",
    )
    return parser.parse_args()


def read_key_values(path):
    values = {}
    with open(path) as handle:
        for line_number, line in enumerate(handle):
            if not line.strip():
                continue
            if "\t" in line:
                key, value = line.rstrip("\n").split("\t", 1)
            else:
                key, value = line.strip().split(None, 1)
            if line_number == 0 and key == "parameter" and value == "value":
                continue
            values[key] = value
    return values


def load_em_parameters(outdir):
    path = os.path.join(outdir, "em_parameters_95.txt")
    if not os.path.exists(path):
        path = os.path.join(outdir, "em_parameters.txt")
    if not os.path.exists(path):
        raise FileNotFoundError("Missing em_parameters_95.txt/em_parameters.txt. Run the reassembly module first.")
    raw = read_key_values(path)
    return {key: float(value) for key, value in raw.items() if key != "parameter"}


def load_pair_counts(outdir, bam, hic1=None, sg1=None):
    path = os.path.join(outdir, "pair_counts_95.txt")
    if not os.path.exists(path):
        path = os.path.join(outdir, "pair_counts.txt")
    counts = {}
    if os.path.exists(path):
        raw = read_key_values(path)
        counts = {key: int(float(value)) for key, value in raw.items()}
    if "all_hic_library_reads" not in counts:
        if not hic1:
            raise ValueError("Missing cached Hi-C count; provide --hic1.")
        counts["all_hic_library_reads"] = reassembly.count_fastq_pairs(hic1)
    if "original_shotgun_library_reads" not in counts:
        if not sg1:
            raise ValueError("Missing cached shotgun count; provide --sg1.")
        counts["original_shotgun_library_reads"] = reassembly.count_fastq_pairs(sg1)
    if "intra_contig_hic_pairs_modeled" not in counts or "inter_contig_hic_pairs" not in counts:
        n, m = reassembly.count_intra_inter_pairs(bam)
        counts["intra_contig_hic_pairs_modeled"] = n
        counts["inter_contig_hic_pairs"] = m
    return counts


def copy_distribution_files(outdir, threshold):
    for stem in ("hic_insert_size", "sg_insert_size"):
        dst = reassembly.suffixed(outdir, stem, threshold)
        if os.path.exists(dst):
            continue
        for src in (
            reassembly.suffixed(outdir, stem, 95),
            os.path.join(outdir, f"{stem}.txt"),
        ):
            if os.path.exists(src):
                shutil.copyfile(src, dst)
                break
    all_hic_dst = reassembly.suffixed(outdir, "all_hic_insert_size", threshold)
    if not os.path.exists(all_hic_dst):
        for src in (
            reassembly.suffixed(outdir, "all_hic_insert_size", 75),
            reassembly.suffixed(outdir, "all_hic_insert_size", 95),
            reassembly.suffixed(outdir, "hic_insert_size", threshold),
            os.path.join(outdir, "all_hic_insert_size.txt"),
            os.path.join(outdir, "hic_insert_size.txt"),
        ):
            if os.path.exists(src):
                shutil.copyfile(src, all_hic_dst)
                break


def run_command(cmd):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def write_threshold_fastqs(args, threshold, readnames):
    suffix = reassembly.threshold_suffix(threshold)
    required = [args.hic1, args.hic2, args.sg1, args.sg2]
    if any(path is None for path in required):
        raise ValueError("--write-fastq requires --hic1, --hic2, --sg1, and --sg2.")

    sg_in_hic_fwd = os.path.join(args.outdir, f"sg_in_hic{suffix}.forward.fastq.gz")
    sg_in_hic_rev = os.path.join(args.outdir, f"sg_in_hic{suffix}.reverse.fastq.gz")
    combined_sg_fwd = os.path.join(args.outdir, f"new_sg_forward{suffix}.fastq.gz")
    combined_sg_rev = os.path.join(args.outdir, f"new_sg_reverse{suffix}.fastq.gz")
    all_readnames = os.path.join(args.outdir, "all_hic_readnames.txt")
    non_sg_readnames = reassembly.suffixed(args.outdir, "readname_non_sg_in_hic", threshold)
    non_sg_fwd = os.path.join(args.outdir, f"non_sg_in_hic{suffix}.forward.fastq.gz")
    non_sg_rev = os.path.join(args.outdir, f"non_sg_in_hic{suffix}.reverse.fastq.gz")

    if not os.path.exists(all_readnames):
        run_command(f"samtools view {args.bam} | cut -f1 | sort -u > {all_readnames}")
    run_command(f"seqtk subseq {args.hic1} {readnames} | gzip > {sg_in_hic_fwd}")
    run_command(f"seqtk subseq {args.hic2} {readnames} | gzip > {sg_in_hic_rev}")
    run_command(f"cat {sg_in_hic_fwd} {args.sg1} > {combined_sg_fwd}")
    run_command(f"cat {sg_in_hic_rev} {args.sg2} > {combined_sg_rev}")
    run_command(f"comm -23 {all_readnames} {readnames} > {non_sg_readnames}")
    run_command(f"seqtk subseq {args.hic1} {non_sg_readnames} | gzip > {non_sg_fwd}")
    run_command(f"seqtk subseq {args.hic2} {non_sg_readnames} | gzip > {non_sg_rev}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    thresholds = [int(item) for item in args.thresholds.split(",") if item.strip()]
    params = load_em_parameters(args.outdir)
    counts = load_pair_counts(args.outdir, args.bam, args.hic1, args.sg1)

    summary = os.path.join(args.outdir, "cutoff_test_summary.tsv")
    with open(summary, "w") as summary_out:
        summary_out.write(
            "threshold\tcutoff\tsg_like_hic_pairs\tactual_hic_fraction\t"
            "estimated_non_informative_fraction\tdifference\n"
        )
        for threshold in thresholds:
            copy_distribution_files(args.outdir, threshold)
            cutoff = reassembly.compute_cutoff(
                params["mu_C"],
                params["mu_N"],
                params["sigma_C"],
                params["sigma_N"],
                threshold,
            )
            cutoff_file = reassembly.suffixed(args.outdir, "hic_insert_size_cutoff", threshold)
            with open(cutoff_file, "w") as out:
                out.write(f"{cutoff:.2f}\n")

            readnames = reassembly.suffixed(args.outdir, "readname_sg_in_hic", threshold)
            sg_like_hic_insert_sizes = reassembly.suffixed(args.outdir, "sg_like_hic_insert_size", threshold)
            sg_like_hic = reassembly.write_sg_like_readnames(
                args.bam,
                cutoff,
                readnames,
                sg_like_hic_insert_sizes,
            )
            fractions = reassembly.suffixed(args.outdir, "fractions", threshold)
            reassembly.write_fraction_report(
                fractions,
                sg_like_hic,
                counts["all_hic_library_reads"],
                counts["original_shotgun_library_reads"],
                counts["intra_contig_hic_pairs_modeled"],
                counts["inter_contig_hic_pairs"],
                params["pi_C"],
                params["pi_N"],
            )
            shutil.copyfile(fractions, reassembly.suffixed(args.outdir, "stats", threshold))
            info = read_key_values(fractions)
            actual = float(info["actual percentage of Hi-C reads used in reassembly"])
            estimated = float(info["estimated non-informative fraction (1 - informative fraction)"])
            summary_out.write(
                f"{threshold}\t{cutoff:.2f}\t{sg_like_hic}\t{actual:.6f}\t"
                f"{estimated:.6f}\t{actual - estimated:.6f}\n"
            )
            if args.write_fastq:
                write_threshold_fastqs(args, threshold, readnames)
    print(f"[INFO] Cutoff test summary written to {summary}")


if __name__ == "__main__":
    main()
