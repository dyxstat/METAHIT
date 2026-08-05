#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import json
import os
import shlex
import subprocess
import sys
from itertools import groupby

import numpy as np
import pysam
from scipy.stats import norm


# ============================================================
# Command-line parsing
# ============================================================

def parse_args():
    script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    default_metahict_path = os.path.abspath(os.path.join(script_dir, "../.."))

    parser = argparse.ArgumentParser(
        description=(
            "EM-based selection of short-insert, WGS-like intra-contig Hi-C "
            "read pairs for reassembly."
        )
    )

    parser.add_argument("--bin", required=True, help="Binning result directory")
    parser.add_argument("--hic1", required=True, help="Hi-C forward FASTQ.gz")
    parser.add_argument("--hic2", required=True, help="Hi-C reverse FASTQ.gz")
    parser.add_argument("--sg1", required=True, help="Shotgun forward FASTQ.gz")
    parser.add_argument("--sg2", required=True, help="Shotgun reverse FASTQ.gz")
    parser.add_argument("--bam", required=True, help="Hi-C BAM mapped to assembly")
    parser.add_argument("--outdir", required=True, help="Output directory")

    parser.add_argument(
        "-p",
        "--metahict_path",
        default=default_metahict_path,
        help="Path to METAHICT folder",
    )

    parser.add_argument("-t", "--threads", type=int, default=80)
    parser.add_argument("-m", "--memory", type=int, default=24)

    parser.add_argument(
        "-k",
        "--top_k",
        type=int,
        default=100,
        help="Number of longest contigs used for EM fitting",
    )

    parser.add_argument(
        "-c",
        "--contig_sort_list",
        required=True,
        help="File listing contig names sorted by length descending",
    )

    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-match-len", type=int, default=30)
    parser.add_argument("--exclude-duplicates", action="store_true")

    parser.add_argument(
        "--cutoff-quantile",
        type=float,
        default=0.95,
        help="Cutoff quantile for N component. Default: 0.95",
    )

    parser.add_argument(
        "--bam-name-sorted",
        action="store_false",
        help=(
            "Set this if --bam is already name-sorted. Otherwise the script "
            "will create a name-sorted BAM in the output directory."
        ),
    )

    parser.add_argument(
        "--write-nonselected-hic",
        action="store_true",
        help="Also write FASTQ files for non-selected Hi-C reads.",
    )

    parser.add_argument(
        "--checkm2_db",
        type=str,
        default=None,
        help="Path to CheckM2 database",
    )
    parser.add_argument("--min-contig-len", type=int, default=500)
    parser.add_argument("--strict-cut-off", type=int, default=2)
    parser.add_argument("--permissive-cut-off", type=int, default=5)
    parser.add_argument("--contamination-penalty", type=float, default=5)
    parser.add_argument("--skip-checkm2", action="store_true")
    parser.add_argument("--tmp-dir", default=None)
    parser.add_argument("--spades-mode", choices=["careful", "none"], default="careful")
    parser.add_argument("--spades-phred-offset", default=None)
    parser.add_argument("--spades-extra-args", default="")
    parser.add_argument("--skip-residual-assembly", action="store_true")
    parser.add_argument("--keep-temp", action="store_true")

    return parser.parse_args()


# ============================================================
# Utilities
# ============================================================

def q(path):
    return shlex.quote(str(path))


def run_command(cmd, check=True):
    print(f"[CMD] {cmd}")
    result = subprocess.run(
        cmd,
        shell=True,
        check=check,
        text=True,
        capture_output=True,
    )
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())
    return result


def open_text(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def count_lines(path):
    n = 0
    with open(path, "r") as f:
        for _ in f:
            n += 1
    return n


def validate_fastq_gz(path, name):
    if not path.endswith(".gz"):
        raise ValueError(f"{name} must be gzip-compressed (.fastq.gz): {path}")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found: {path}")


# ============================================================
# BAM pair filtering and mapped insert-size extraction
# ============================================================

def matched_length(aln):
    """
    Count aligned matching operations from CIGAR.

    CIGAR:
      0 = M
      7 = =
      8 = X
    """
    if aln.cigartuples is None:
        return 0
    return sum(length for op, length in aln.cigartuples if op in (0, 7, 8))


def passes_filter(aln, min_mapq, min_match_len, exclude_duplicates):
    if aln.is_unmapped:
        return False
    if aln.is_secondary or aln.is_supplementary:
        return False
    if exclude_duplicates and aln.is_duplicate:
        return False
    if aln.mapping_quality < min_mapq:
        return False
    if aln.reference_id < 0:
        return False
    if aln.reference_start is None or aln.reference_end is None:
        return False
    if matched_length(aln) < min_match_len:
        return False
    return True


def choose_best(alns):
    return sorted(
        alns,
        key=lambda a: (
            a.mapping_quality,
            matched_length(a),
            a.query_alignment_length or 0,
        ),
        reverse=True,
    )[0]


def process_group(reads, min_mapq, min_match_len, exclude_duplicates):
    """
    Return one filtered read1/read2 pair from a query group.
    """
    passing = [
        r for r in reads
        if passes_filter(
            r,
            min_mapq=min_mapq,
            min_match_len=min_match_len,
            exclude_duplicates=exclude_duplicates,
        )
    ]

    if len(passing) < 2:
        return None, "no_two_passing_mates"

    r1 = [r for r in passing if r.is_read1]
    r2 = [r for r in passing if r.is_read2]

    if r1 and r2:
        return (choose_best(r1), choose_best(r2)), "paired"

    # Fallback for unusual BAMs missing read1/read2 flags
    best = sorted(
        passing,
        key=lambda a: (
            a.mapping_quality,
            matched_length(a),
            a.query_alignment_length or 0,
        ),
        reverse=True,
    )[:2]

    if len(best) == 2:
        return (best[0], best[1]), "fallback_two_alignments"

    return None, "no_pair"


def mapped_insert_size(a, b):
    """
    Mapped insert size / mapped template span.

    This is NOT SAM TLEN and NOT the inner gap.

    Definition:
        d = max(end1, end2) - min(start1, start2)
    """
    if a.reference_id != b.reference_id:
        return None

    if a.reference_start is None or a.reference_end is None:
        return None
    if b.reference_start is None or b.reference_end is None:
        return None

    d = max(a.reference_end, b.reference_end) - min(
        a.reference_start,
        b.reference_start,
    )

    if d <= 0:
        return None

    return int(d)


def load_top_contigs(contig_sort_list, top_k):
    if top_k <= 0:
        raise ValueError("--top_k must be > 0")

    contigs = []

    with open(contig_sort_list, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # First column is contig name
            contigs.append(line.split()[0])

            if len(contigs) >= top_k:
                break

    if len(contigs) == 0:
        raise ValueError("No contigs loaded from --contig_sort_list")

    return contigs, set(contigs)


def ensure_name_sorted_bam(bam_path, outdir, threads, already_name_sorted=False):
    if already_name_sorted:
        return bam_path

    name_sorted_bam = os.path.join(outdir, "hic.name_sorted.bam")

    if os.path.exists(name_sorted_bam) and os.path.getsize(name_sorted_bam) > 0:
        print(f"[INFO] Using existing name-sorted BAM: {name_sorted_bam}")
        return name_sorted_bam

    print("[INFO] Creating name-sorted BAM for pair-level extraction")
    run_command(
        f"samtools sort -n -@ {threads} -o {q(name_sorted_bam)} {q(bam_path)}"
    )

    return name_sorted_bam


def extract_insert_size_tables(
    bam_path,
    all_out,
    top_out,
    counts_out,
    top_contig_set,
    min_mapq=30,
    min_match_len=30,
    exclude_duplicates=False,
):
    """
    Traverse a name-sorted BAM once and write:

      1. all intra-contig mapped insert sizes
      2. top-contig intra-contig mapped insert sizes for EM fitting

    Both files have:
      read_id, contig, d, start1, end1, start2, end2, mapq1, mapq2, match_len1, match_len2
    """

    ensure_parent_dir(all_out)
    ensure_parent_dir(top_out)
    ensure_parent_dir(counts_out)

    counts = {
        "bam": bam_path,
        "query_groups": 0,
        "filtered_pairs": 0,
        "intra_pairs": 0,
        "inter_pairs": 0,
        "all_intra_pairs_written": 0,
        "top_intra_pairs_written": 0,
        "top_intra_pairs_skipped": 0,
        "invalid_insert_size_pairs": 0,
        "no_two_passing_mates": 0,
        "fallback_two_alignments": 0,
        "min_mapq": min_mapq,
        "min_match_len": min_match_len,
        "exclude_duplicates": exclude_duplicates,
        "distance_definition": "mapped_insert_size=max(end1,end2)-min(start1,start2)",
        "top_contig_count": len(top_contig_set),
    }

    header = (
        "read_id\tcontig\td\tstart1\tend1\tstart2\tend2\t"
        "mapq1\tmapq2\tmatch_len1\tmatch_len2\n"
    )

    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        open_text(all_out, "wt") as all_f,
        open_text(top_out, "wt") as top_f,
    ):
        all_f.write(header)
        top_f.write(header)

        for qname, group_iter in groupby(
            bam.fetch(until_eof=True),
            key=lambda x: x.query_name,
        ):
            reads = list(group_iter)
            counts["query_groups"] += 1

            pair, status = process_group(
                reads,
                min_mapq=min_mapq,
                min_match_len=min_match_len,
                exclude_duplicates=exclude_duplicates,
            )

            if pair is None:
                counts["no_two_passing_mates"] += 1
                continue

            if status == "fallback_two_alignments":
                counts["fallback_two_alignments"] += 1

            a, b = pair
            counts["filtered_pairs"] += 1

            if a.reference_id != b.reference_id:
                counts["inter_pairs"] += 1
                continue

            counts["intra_pairs"] += 1

            d = mapped_insert_size(a, b)

            if d is None:
                counts["invalid_insert_size_pairs"] += 1
                continue

            contig = bam.get_reference_name(a.reference_id)

            row = (
                f"{qname}\t{contig}\t{d}\t"
                f"{a.reference_start}\t{a.reference_end}\t"
                f"{b.reference_start}\t{b.reference_end}\t"
                f"{a.mapping_quality}\t{b.mapping_quality}\t"
                f"{matched_length(a)}\t{matched_length(b)}\n"
            )

            all_f.write(row)
            counts["all_intra_pairs_written"] += 1

            if contig in top_contig_set:
                top_f.write(row)
                counts["top_intra_pairs_written"] += 1
            else:
                counts["top_intra_pairs_skipped"] += 1

    with open(counts_out, "w") as f:
        json.dump(counts, f, indent=2)

    print(json.dumps(counts, indent=2))

    return counts


# ============================================================
# EM model
# ============================================================

def load_insert_sizes_from_table(path, col="d"):
    values = []

    with open_text(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")

        if col not in header:
            raise ValueError(f"Column '{col}' not found in {path}")

        d_idx = header.index(col)

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            try:
                d = float(fields[d_idx])
            except Exception:
                continue

            if np.isfinite(d) and d > 0:
                values.append(d)

    return np.asarray(values, dtype=float)


def init_params(data, init_frac=0.8):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data) & (data > 0)]

    if len(data) == 0:
        raise ValueError("No valid insert-size data for EM fitting.")

    qv = np.percentile(data, init_frac * 100)

    lower = data[data <= qv]
    upper = data[data > qv]

    if len(lower) == 0:
        lower = data[:max(1, len(data) // 2)]
    if len(upper) == 0:
        upper = data[max(1, len(data) // 2):]

    mu_N = lower.mean()
    mu_C = upper.mean()

    sigma_N = max(lower.std(ddof=1), 1.0) if len(lower) > 1 else 1.0
    sigma_C = max(upper.std(ddof=1), 1.0) if len(upper) > 1 else 1.0

    pi_N = len(lower) / len(data)
    pi_C = 1.0 - pi_N

    return mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C


def em_mix(data, mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C, tol=1e-2, max_iter=100):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data) & (data > 0)]

    if len(data) == 0:
        raise ValueError("No data for EM algorithm.")

    logL = None

    for iteration in range(max_iter):
        wC = pi_C * norm.pdf(data, mu_C, sigma_C)
        wN = pi_N * norm.pdf(data, mu_N, sigma_N)

        total = np.maximum(wC + wN, 1e-10)

        gamma_C = wC / total
        gamma_N = wN / total

        sum_C = gamma_C.sum()
        sum_N = gamma_N.sum()

        if sum_C > 0:
            mu_C = np.sum(gamma_C * data) / sum_C
            sigma_C = max(
                np.sqrt(np.sum(gamma_C * (data - mu_C) ** 2) / sum_C),
                1e-6,
            )

        if sum_N > 0:
            mu_N = np.sum(gamma_N * data) / sum_N
            sigma_N = max(
                np.sqrt(np.sum(gamma_N * (data - mu_N) ** 2) / sum_N),
                1e-6,
            )

        pi_C = gamma_C.mean()
        pi_N = gamma_N.mean()

        density = (
            pi_C * norm.pdf(data, mu_C, sigma_C)
            + pi_N * norm.pdf(data, mu_N, sigma_N)
        )

        new_logL = np.sum(np.log(np.maximum(density, 1e-300)))

        if logL is not None and abs(new_logL - logL) < tol:
            print(f"[INFO] EM converged after {iteration + 1} iterations")
            break

        logL = new_logL

        print(
            f"Iter {iteration + 1}: "
            f"logL={new_logL:.2f}, "
            f"pi_C={pi_C:.4f}, pi_N={pi_N:.4f}, "
            f"mu_C={mu_C:.2f}, mu_N={mu_N:.2f}"
        )

    return mu_C, mu_N, sigma_C, sigma_N, pi_C, pi_N


def fit_em(data, cutoff_quantile=0.95):
    mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C = init_params(data)

    mu_C, mu_N, sigma_C, sigma_N, pi_C, pi_N = em_mix(
        data,
        mu_N,
        mu_C,
        sigma_N,
        sigma_C,
        pi_N,
        pi_C,
    )

    cutoff = mu_N + norm.ppf(cutoff_quantile) * sigma_N

    return {
        "mu_N": float(mu_N),
        "sigma_N": float(sigma_N),
        "pi_N": float(pi_N),
        "mu_C": float(mu_C),
        "sigma_C": float(sigma_C),
        "pi_C": float(pi_C),
        "cutoff": float(cutoff),
        "cutoff_quantile": float(cutoff_quantile),
    }


# ============================================================
# Read-name selection
# ============================================================

def write_selected_readnames_from_insert_table(insert_table, cutoff, out_readnames):
    """
    Select read names from the same pair-level insert-size table used by the method.

    This avoids SAM TLEN and avoids selecting reads from alignments that were not
    part of the filtered intra-contig pair set.
    """

    ensure_parent_dir(out_readnames)

    tmp = out_readnames + ".tmp"

    selected = 0
    total = 0

    with open_text(insert_table, "rt") as f, open(tmp, "w") as out:
        header = f.readline().rstrip("\n").split("\t")

        if "read_id" not in header:
            raise ValueError("read_id column missing from insert-size table.")
        if "d" not in header:
            raise ValueError("d column missing from insert-size table.")

        rid_idx = header.index("read_id")
        d_idx = header.index("d")

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            try:
                d = float(fields[d_idx])
            except Exception:
                continue

            if not np.isfinite(d) or d <= 0:
                continue

            total += 1

            if d <= cutoff:
                out.write(fields[rid_idx] + "\n")
                selected += 1

    run_command(f"LC_ALL=C sort -u {q(tmp)} > {q(out_readnames)}")
    os.remove(tmp)

    selected_unique = count_lines(out_readnames)

    return selected_unique, total


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.checkm2_db:
        os.environ["CHECKM2DB"] = args.checkm2_db
        print(f"[INFO] Using CheckM2 database: {args.checkm2_db}")
    else:
        print("[INFO] Using default CheckM2 database path.")

    validate_fastq_gz(args.hic1, "--hic1")
    validate_fastq_gz(args.hic2, "--hic2")
    validate_fastq_gz(args.sg1, "--sg1")
    validate_fastq_gz(args.sg2, "--sg2")

    # ------------------------------------------------------------
    # Prepare name-sorted BAM
    # ------------------------------------------------------------
    name_sorted_bam = ensure_name_sorted_bam(
        bam_path=args.bam,
        outdir=args.outdir,
        threads=args.threads,
        already_name_sorted=args.bam_name_sorted,
    )

    # ------------------------------------------------------------
    # Load top contigs
    # ------------------------------------------------------------
    top_contigs, top_contig_set = load_top_contigs(
        args.contig_sort_list,
        args.top_k,
    )

    top_contigs_file = os.path.join(args.outdir, f"top_{args.top_k}_contigs.txt")
    with open(top_contigs_file, "w") as f:
        for c in top_contigs:
            f.write(c + "\n")

    print(f"[INFO] Top {args.top_k} contigs written to {top_contigs_file}")

    # ------------------------------------------------------------
    # Extract mapped insert sizes
    # ------------------------------------------------------------
    all_insert_table = os.path.join(args.outdir, "all_intra_insert_sizes.tsv.gz")
    top_insert_table = os.path.join(args.outdir, f"top{args.top_k}_intra_insert_sizes.tsv.gz")
    counts_json = os.path.join(args.outdir, "insert_size_counts.json")

    counts = extract_insert_size_tables(
        bam_path=name_sorted_bam,
        all_out=all_insert_table,
        top_out=top_insert_table,
        counts_out=counts_json,
        top_contig_set=top_contig_set,
        min_mapq=args.min_mapq,
        min_match_len=args.min_match_len,
        exclude_duplicates=args.exclude_duplicates,
    )

    # ------------------------------------------------------------
    # Fit EM on top-k insert sizes
    # ------------------------------------------------------------
    print(f"[INFO] Fitting EM using insert sizes from top {args.top_k} contigs")

    top_insert_sizes = load_insert_sizes_from_table(top_insert_table, col="d")

    if len(top_insert_sizes) == 0:
        raise RuntimeError("No top-contig insert sizes available for EM fitting.")

    em_params = fit_em(
        top_insert_sizes,
        cutoff_quantile=args.cutoff_quantile,
    )

    cutoff = em_params["cutoff"]

    print("[INFO] EM parameters:")
    print(json.dumps(em_params, indent=2))
    print(f"[INFO] Insert-size cutoff: {cutoff:.2f}")

    with open(os.path.join(args.outdir, "em_parameters.json"), "w") as f:
        json.dump(em_params, f, indent=2)

    with open(os.path.join(args.outdir, "mixing_proportion.txt"), "w") as f:
        f.write(f"{em_params['pi_C']}\n")

    with open(os.path.join(args.outdir, "insert_size_cutoff.txt"), "w") as f:
        f.write(f"{cutoff:.2f}\n")

    # ------------------------------------------------------------
    # Library-level diagnostics
    # ------------------------------------------------------------
    n = counts["intra_pairs"]
    m = counts["inter_pairs"]
    total_pair_count = n + m

    informative_fraction = (
        (m + n * em_params["pi_C"]) / total_pair_count
        if total_pair_count > 0
        else float("nan")
    )

    ratio_3d = (
        m / total_pair_count
        if total_pair_count > 0
        else float("nan")
    )

    with open(os.path.join(args.outdir, "informative_fraction.txt"), "w") as f:
        f.write(f"{informative_fraction}\n")

    # Backward-compatible filename
    with open(os.path.join(args.outdir, "long_range_ratio.txt"), "w") as f:
        f.write(f"{informative_fraction}\n")

    with open(os.path.join(args.outdir, "3d_ratio.txt"), "w") as f:
        f.write(f"{ratio_3d:.6f}\n")

    # ------------------------------------------------------------
    # Select rescued Hi-C read names from all intra-contig pairs
    # ------------------------------------------------------------
    readnames = os.path.join(args.outdir, "readname_sg_in_hic.txt")

    selected_pairs, all_intra_written = write_selected_readnames_from_insert_table(
        insert_table=all_insert_table,
        cutoff=cutoff,
        out_readnames=readnames,
    )

    selected_fraction_intra = (
        selected_pairs / all_intra_written
        if all_intra_written > 0
        else float("nan")
    )

    selected_fraction_filtered = (
        selected_pairs / counts["filtered_pairs"]
        if counts["filtered_pairs"] > 0
        else float("nan")
    )

    summary = {
        "top_k": args.top_k,
        "top_insert_sizes_for_EM": int(len(top_insert_sizes)),
        "all_intra_insert_pairs": int(all_intra_written),
        "selected_pairs": int(selected_pairs),
        "selected_fraction_of_intra_pairs": selected_fraction_intra,
        "selected_fraction_of_filtered_pairs": selected_fraction_filtered,
        "informative_fraction": informative_fraction,
        "one_minus_informative_fraction": (
            1 - informative_fraction
            if np.isfinite(informative_fraction)
            else float("nan")
        ),
        "3d_ratio": ratio_3d,
        **em_params,
    }

    with open(os.path.join(args.outdir, "em_selection_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("[INFO] EM selection summary:")
    print(json.dumps(summary, indent=2))

    # ------------------------------------------------------------
    # Extract selected Hi-C reads and merge with shotgun reads
    # ------------------------------------------------------------
    sg_in_hic_fwd = os.path.join(args.outdir, "sg_in_hic.forward.fastq.gz")
    sg_in_hic_rev = os.path.join(args.outdir, "sg_in_hic.reverse.fastq.gz")

    combined_sg_fwd = os.path.join(args.outdir, "new_sg_forward.fastq.gz")
    combined_sg_rev = os.path.join(args.outdir, "new_sg_reverse.fastq.gz")

    run_command(
        f"seqtk subseq {q(args.hic1)} {q(readnames)} | gzip -c > {q(sg_in_hic_fwd)}"
    )
    run_command(
        f"seqtk subseq {q(args.hic2)} {q(readnames)} | gzip -c > {q(sg_in_hic_rev)}"
    )

    run_command(f"cat {q(sg_in_hic_fwd)} {q(args.sg1)} > {q(combined_sg_fwd)}")
    run_command(f"cat {q(sg_in_hic_rev)} {q(args.sg2)} > {q(combined_sg_rev)}")

    # ------------------------------------------------------------
    # Optional: output non-selected Hi-C reads
    # ------------------------------------------------------------
    if args.write_nonselected_hic:
        all_readnames = os.path.join(args.outdir, "all_hic_readnames.txt")
        non_selected_readnames = os.path.join(args.outdir, "readname_non_sg_in_hic.txt")

        run_command(
            f"samtools view -F 0x900 {q(name_sorted_bam)} "
            f"| cut -f1 | LC_ALL=C sort -u > {q(all_readnames)}"
        )

        run_command(
            f"LC_ALL=C comm -23 {q(all_readnames)} {q(readnames)} "
            f"> {q(non_selected_readnames)}"
        )

        non_sg_fwd = os.path.join(args.outdir, "non_sg_in_hic.forward.fastq.gz")
        non_sg_rev = os.path.join(args.outdir, "non_sg_in_hic.reverse.fastq.gz")

        run_command(
            f"seqtk subseq {q(args.hic1)} {q(non_selected_readnames)} | gzip -c > {q(non_sg_fwd)}"
        )
        run_command(
            f"seqtk subseq {q(args.hic2)} {q(non_selected_readnames)} | gzip -c > {q(non_sg_rev)}"
        )

    # ------------------------------------------------------------
    # Reassemble bins
    # ------------------------------------------------------------
    reassemble_cmd = (
        f"bash {q(os.path.join(args.metahict_path, '7_reassembly/scripts/reassemble_bins.sh'))} "
        f"-b {q(args.bin)} "
        f"-o {q(args.outdir)} "
        f"-1 {q(combined_sg_fwd)} "
        f"-2 {q(combined_sg_rev)} "
        f"-m {args.memory} "
        f"-t {args.threads} "
        f"-l {args.min_contig_len} "
        f"--strict-cut-off {args.strict_cut_off} "
        f"--permissive-cut-off {args.permissive_cut_off} "
        f"--contamination-penalty {args.contamination_penalty} "
        f"--spades-mode {q(args.spades_mode)} "
        f"--parallel"
    )

    if args.checkm2_db:
        reassemble_cmd += f" --checkm2_db {q(args.checkm2_db)}"
    if args.skip_checkm2:
        reassemble_cmd += " --skip-checkm"
    if args.tmp_dir:
        reassemble_cmd += f" --tmp-dir {q(args.tmp_dir)}"
    if args.spades_phred_offset:
        reassemble_cmd += f" --spades-phred-offset {q(args.spades_phred_offset)}"
    if args.spades_extra_args:
        reassemble_cmd += f" --spades-extra-args={q(args.spades_extra_args)}"
    if args.skip_residual_assembly:
        reassemble_cmd += " --skip-residual-assembly"
    if args.keep_temp:
        reassemble_cmd += " --keep-temp"

    run_command(reassemble_cmd)

    print("[INFO] Reassembly pipeline completed successfully!")


if __name__ == "__main__":
    main()
