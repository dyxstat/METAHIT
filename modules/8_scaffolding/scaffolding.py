#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scaffolding Pipeline for Metagenomic Genome Assembly

Implements full scaffolding workflow:
1. Filter input FASTA (<5kb removed)
2. FASTA indexing
3. Align Hi-C reads or use provided BAM
4. Run YaHS scaffolding
5. Cut scaffolds into segments
6. Realign Hi-C to scaffolded genome
7. Generate MetaCC contact matrix
8. Normalize with NormCC
9. Create scaffold mapping
10. Draw heatmap
11. Compute metrics (N50/L50, CheckM2, enrichment)
"""

import argparse
import logging
import sys
import os
import subprocess
import shlex
import shutil
import warnings
import gzip
import pickle
import pandas as pd
import scipy.sparse as scisp
from pandas.errors import SettingWithCopyWarning
warnings.filterwarnings("ignore", category=SettingWithCopyWarning)

from MetaCC.Script.normalized_contact import NormCCMap
from MetaCC.Script.utils import save_object, make_dir
from MetaCC.Script.normcc import normcc
from scaffolding_eligibility import (
    assess_scaffolding_eligibility,
    fasta_sequence_lengths,
    write_scaffolding_status,
)

class ApplicationException(Exception):
    pass


def filter_fasta_by_length(input_fasta, output_fasta, min_length=5000):
    retained = 0
    seq, header = "", ""
    with open(input_fasta) as f_in, open(output_fasta, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                if header and len(seq) >= min_length:
                    f_out.write(header + "\n" + seq + "\n")
                    retained += 1
                header, seq = line.strip(), ""
            else:
                seq += line.strip()
        if header and len(seq) >= min_length:
            f_out.write(header + "\n" + seq + "\n")
            retained += 1
    return retained


def validate_bam_reference(bam, fasta):
    """Require every scaffolding FASTA sequence to match the BAM reference."""
    import pysam

    fasta_lengths = fasta_sequence_lengths(fasta)
    try:
        with pysam.AlignmentFile(bam, "rb") as handle:
            bam_lengths = dict(zip(handle.references, handle.lengths))
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot read the supplied scaffolding BAM {bam}: {error}") from error

    missing = sorted(set(fasta_lengths) - set(bam_lengths))
    mismatched = sorted(
        name
        for name, length in fasta_lengths.items()
        if name in bam_lengths and bam_lengths[name] != length
    )
    if missing or mismatched:
        details = []
        if missing:
            details.append("missing reference(s): " + ", ".join(missing[:10]))
        if mismatched:
            details.append("length mismatch(es): " + ", ".join(mismatched[:10]))
        raise ValueError(
            "The supplied --bam was not aligned to the selected scaffolding bin; "
            + "; ".join(details)
        )

    extra = sorted(set(bam_lengths) - set(fasta_lengths))
    if extra:
        logger.warning(
            "The supplied BAM contains %d reference(s) excluded by --min-contig-len; "
            "only references retained in the filtered bin will contribute",
            len(extra),
        )


def filter_metacc_inputs(contact_matrix, contig_info):
    """Apply the original MetaCC site/coverage eligibility filters consistently."""
    required = {"name", "sites", "covcc"}
    missing = sorted(required - set(contig_info.columns))
    if missing:
        raise ValueError(
            "MetaCC contig table is missing required column(s): " + ", ".join(missing)
        )
    if contact_matrix.shape[0] != len(contig_info):
        raise ValueError(
            f"Raw MetaCC matrix dimension {contact_matrix.shape[0]} does not match "
            f"the {len(contig_info)} contig rows"
        )
    keep = (contig_info["sites"] > 0).to_numpy() & (contig_info["covcc"] > 0).to_numpy()
    if not keep.any():
        raise ValueError("No MetaCC contigs remain after restriction-site and coverage filtering")
    filtered_matrix = contact_matrix.tocsr()[keep, :][:, keep]
    filtered_info = contig_info.loc[keep].reset_index(drop=True)
    signal_matrix = filtered_matrix.copy()
    signal_matrix.setdiag(0)
    signal_matrix.eliminate_zeros()
    filtered_info["signal"] = signal_matrix.sum(axis=1).A1
    return filtered_matrix, filtered_info


def cut_scaffolds_to_segments(fasta_in, fasta_out, seg_len):
    mapping = {}
    with open(fasta_in) as f_in, open(fasta_out, "w") as f_out:
        scaffold, seq = None, ""
        for line in f_in:
            if line.startswith(">"):
                if scaffold and seq:
                    for i in range(0, len(seq), seg_len):
                        seg = seq[i:i+seg_len]
                        if seg:
                            seg_name = f"{scaffold}_seg_{i//seg_len+1}"
                            f_out.write(f">{seg_name}\n{seg}\n")
                            mapping[seg_name] = scaffold
                scaffold, seq = line[1:].strip().split()[0], ""
            else:
                seq += line.strip()
        if scaffold and seq:
            for i in range(0, len(seq), seg_len):
                seg = seq[i:i+seg_len]
                if seg:
                    seg_name = f"{scaffold}_seg_{i//seg_len+1}"
                    f_out.write(f">{seg_name}\n{seg}\n")
                    mapping[seg_name] = scaffold
    return mapping


def create_scaffold_mapping(mapping, outpath):
    save_object(outpath.replace(".gz", ""), mapping)
    with gzip.open(outpath, "wb") as f:
        pickle.dump(mapping, f)
    return mapping


def calculate_n50_l50(fasta):
    lens, cur = [], 0
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lens.append(cur)
    if not lens:
        return 0, 0
    lens.sort(reverse=True)
    total, half, cum = sum(lens), sum(lens)/2, 0
    for i, l in enumerate(lens):
        cum += l
        if cum >= half:
            return l, i+1
    return 0, 0


def run_checkm2(fasta, outdir, threads=1, checkm2_db=None):
    env = os.environ.copy()
    if checkm2_db:
        env["CHECKM2DB"] = checkm2_db
    executable = os.environ.get("METAHICT_CHECKM2_BIN") or shutil.which("checkm2")
    if not executable:
        logger.warning("CheckM2 was not found; run ./metahict install or set METAHICT_CHECKM2_BIN")
        return None, None
    executable_dir = os.path.dirname(os.path.abspath(os.path.expanduser(executable)))
    path_entries = [
        entry
        for entry in env.get("PATH", "").split(os.pathsep)
        if entry and entry != executable_dir
    ]
    env["PATH"] = os.pathsep.join([executable_dir, *path_entries])
    if not shutil.which("diamond", path=env["PATH"]):
        logger.warning("DIAMOND was not found in the CheckM2 environment: %s", executable_dir)
        return None, None
    cmd = [
        executable,
        "predict",
        "--threads", str(threads),
        "--input", str(fasta),
        "--output-directory", str(outdir),
    ]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if res.returncode != 0:
            logger.warning(f"CheckM2 failed: {res.stderr}")
            return None, None
        qfile = os.path.join(outdir, "quality_report.tsv")
        if os.path.exists(qfile):
            df = pd.read_csv(qfile, sep="\t")
            if not df.empty:
                return df["Completeness"].iloc[0], df["Contamination"].iloc[0]
    except Exception as e:
        logger.warning(f"CheckM2 error: {e}")
    return None, None


def calculate_hic_enrichment_contigs(bam, fasta):
    import pysam
    contigs = {l[1:].strip().split()[0] for l in open(fasta) if l.startswith(">")}
    within, across = 0, 0
    with pysam.AlignmentFile(bam, "rb") as bf:
        for r in bf.fetch(until_eof=True):
            if r.is_paired and not r.is_unmapped and not r.mate_is_unmapped:
                c1, c2 = r.reference_name, r.next_reference_name
                if c1 in contigs and c2 in contigs:
                    if c1 == c2: within += 1
                    else: across += 1
    return within / across if across else (float("inf") if within else 0)


def calculate_hic_enrichment_scaffolds(bam, fasta):
    import pysam
    scaffs = {l[1:].strip().split()[0] for l in open(fasta) if l.startswith(">")}
    def sname(seg): return seg.split("_seg_")[0] if "_seg_" in seg else seg
    within, across = 0, 0
    with pysam.AlignmentFile(bam, "rb") as bf:
        for r in bf.fetch(until_eof=True):
            if r.is_paired and not r.is_unmapped and not r.mate_is_unmapped:
                s1, s2 = sname(r.reference_name), sname(r.next_reference_name)
                if s1 in scaffs and s2 in scaffs:
                    if s1 == s2: within += 1
                    else: across += 1
    return within / across if across else (float("inf") if within else 0)


def calculating_metrics(filt_fa, genome, args, original_bam, scaffold_bam, quality_dir):
    logger.info("Calculating quality metrics")
    n50_o, l50_o = calculate_n50_l50(filt_fa)
    n50_s, l50_s = calculate_n50_l50(genome)
    c_o = sum(1 for l in open(filt_fa) if l.startswith(">"))
    c_s = sum(1 for l in open(genome) if l.startswith(">"))
    t = int(args.t) if args.t else 1
    comp_o = cont_o = comp_s = cont_s = None
    if not args.skip_checkm2:
        d1 = os.path.join(quality_dir, "checkm2_original")
        d2 = os.path.join(quality_dir, "checkm2_scaffolded")
        make_dir(d1); make_dir(d2)
        db = getattr(args, "checkm2_db", None)
        comp_o, cont_o = run_checkm2(filt_fa, d1, t, db)
        comp_s, cont_s = run_checkm2(genome, d2, t, db)
    enr_o, enr_s = None, None
    if os.path.exists(original_bam):
        enr_o = calculate_hic_enrichment_contigs(original_bam, filt_fa)
    if os.path.exists(scaffold_bam):
        enr_s = calculate_hic_enrichment_scaffolds(scaffold_bam, genome)
    return {
        "Original_N50": n50_o, "Original_L50": l50_o, "Original_contig_count": c_o,
        "Scaffolded_N50": n50_s, "Scaffolded_L50": l50_s, "Scaffolded_scaffold_count": c_s,
        "Original_completeness": comp_o, "Original_contamination": cont_o,
        "Scaffolded_completeness": comp_s, "Scaffolded_contamination": cont_s,
        "Original_HiC_enrichment_ratio": enr_o, "Scaffolded_HiC_enrichment_ratio": enr_s
    }


def save_metrics(metrics, outfile):
    with open(outfile, "w") as f:
        f.write("Scaffolding Quality Metrics\n==========================\n\n")
        for k, v in metrics.items():
            f.write(f"{k}: {v if v is not None else 'N/A'}\n")


def cleanup_intermediates(intermediate_dir):
    if os.path.isdir(intermediate_dir):
        shutil.rmtree(intermediate_dir)


def run_command(argv):
    logger.info("Executing: %s", shlex.join([str(value) for value in argv]))
    subprocess.run([str(value) for value in argv], check=True)


def alignment_command(args, reference, output_dir):
    """Build the shared Module 3 Python alignment command.

    Values that begin with ``-`` must use argparse's ``--option=value`` form.
    Passing them as a following argv item makes argparse interpret values such
    as ``-5SP`` and ``-F 0x900`` as new options instead of option values.
    """
    return [
        sys.executable,
        os.path.join(args.p, "3_alignment", "alignment.py"),
        "-p", args.p,
        "-r", reference,
        "-1", args.hic1,
        "-2", args.hic2,
        "-o", output_dir,
        "-t", str(args.t),
        f"--bwa-options={args.bwa_options}",
        f"--samtools-filter={args.samtools_filter}",
        "--sort-memory", args.sort_memory,
        "--mapq", "0",
        "--min-match-len", "0",
        "--skip-metrics",
    ]


def parse_restriction_enzymes(value):
    """Return a validated, whitespace-normalized restriction-enzyme list."""
    enzymes = [enzyme.strip() for enzyme in value.split(",") if enzyme.strip()]
    if not enzymes:
        raise argparse.ArgumentTypeError(
            "at least one restriction enzyme is required, for example Sau3AI,MluCI"
        )
    return enzymes


def build_parser():
    parser = argparse.ArgumentParser(
        description="Scaffold one MAG with Hi-C contacts and YaHS",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Routine use: ./metahict run --entry-module scaffolding --help",
    )
    parser.add_argument("-p", required=True, help="Base METAHICT modules path")
    parser.add_argument("--fasta", required=True, help="FASTA for the one bin to scaffold")
    parser.add_argument(
        "--enzyme",
        required=True,
        type=parse_restriction_enzymes,
        help="Comma-separated restriction enzyme list",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--hic1", required=True, help="Cleaned Hi-C forward FASTQ")
    parser.add_argument("--hic2", required=True, help="Cleaned Hi-C reverse FASTQ")

    parser.add_argument("--bam", help="Optional Hi-C BAM aligned to the selected bin")
    parser.add_argument(
        "-t", "--threads", dest="t", type=int, default=8,
        help="Threads passed to internal alignment, SAMtools, and CheckM2",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        dest="r",
        type=int,
        default=10000,
        help="Heatmap segment resolution",
    )
    parser.add_argument("--min-contig-len", type=int, default=5000, help="Minimum contig length retained for scaffolding")
    parser.add_argument("--bwa-options", default="-5SP", help="BWA MEM options for internal alignments")
    parser.add_argument("--samtools-filter", default="-F 0x900", help="samtools view filter for internal alignments")
    parser.add_argument("--sort-memory", default="1G", help="Memory per samtools sort thread")
    parser.add_argument("--metacc-min-mapq", type=int, default=30, help="Minimum MAPQ for MetaCC contact generation")
    parser.add_argument("--metacc-min-len", type=int, default=1000, help="Minimum contig length for MetaCC contact generation")
    parser.add_argument("--metacc-min-match", type=int, default=30, help="Minimum aligned match length for MetaCC contact generation")
    parser.add_argument("--metacc-min-signal", type=int, default=2, help="Minimum signal for MetaCC contact generation")
    parser.add_argument("--yahs-resolutions", default="", help="YaHS scaffolding resolution list; empty uses YaHS automatic selection")
    parser.add_argument("--yahs-min-mapq", type=int, default=10, help="YaHS minimum mapping quality")
    parser.add_argument("--yahs-min-contig-len", type=int, default=0, help="YaHS minimum contig length to scaffold")
    parser.add_argument("--yahs-rounds", type=int, default=1, help="YaHS rounds at each resolution level")
    parser.add_argument("--yahs-no-contig-ec", action="store_true", help="Disable YaHS contig error correction")
    parser.add_argument("--yahs-no-scaffold-ec", action="store_true", help="Disable YaHS scaffold error correction")
    parser.add_argument("--yahs-no-mem-check", action="store_true", help="Disable YaHS runtime memory check")
    parser.add_argument("--yahs-extra-args", default="", help="Additional options passed directly to YaHS")
    parser.add_argument("--normcc-thres", type=float, default=0.05, help="NormCC denoising threshold")
    parser.add_argument("--heatmap-max-image", type=int, default=5000, help="Maximum heatmap image dimension before downsampling")
    parser.add_argument("--skip-checkm2", action="store_true", help="Skip CheckM2 quality evaluation")
    parser.add_argument("--checkm2_db", help="Path to custom CheckM2 database")
    parser.add_argument("--tmp-dir", default=os.environ.get("METAHICT_TMP_ROOT", os.environ.get("TMPDIR", "/tmp")),
                        help="Temporary directory root")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files for debugging")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()

    if args.t < 1:
        parser.error("--threads must be at least 1")
    if args.r < 1 or args.heatmap_max_image < 1 or args.yahs_rounds < 1:
        parser.error("heatmap resolution, heatmap size, and YaHS rounds must be positive")
    if min(
        args.min_contig_len,
        args.metacc_min_mapq,
        args.metacc_min_len,
        args.metacc_min_match,
        args.metacc_min_signal,
        args.yahs_min_mapq,
        args.yahs_min_contig_len,
    ) < 0:
        parser.error("scaffolding and contact filtering thresholds must be non-negative")
    if not 0 <= args.normcc_thres <= 1:
        parser.error("--normcc-thres must be between 0 and 1")

    if os.path.isdir(os.path.join(args.p, "modules", "8_scaffolding")):
        args.p = os.path.join(args.p, "modules")
    if not os.path.isfile(os.path.join(args.p, "8_scaffolding", "heatmap.py")):
        parser.error(f"cannot locate the flat 8_scaffolding module from -p {args.p}")

    make_dir(args.outdir)
    make_dir(args.tmp_dir)
    os.environ["TMPDIR"] = args.tmp_dir
    log_file = os.path.join(args.outdir, "scaffolding.log")
    logging.basicConfig(level=logging.INFO, format="%(levelname)-8s | %(asctime)s | %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="a")])
    global logger
    logger = logging.getLogger("main")

    try:
        intermediate_dir = os.path.join(args.outdir, "intermediates")
        make_dir(intermediate_dir)
        # Step 1: filter
        filt_fa = os.path.join(intermediate_dir, "filtered_bin.fa")
        assessment = assess_scaffolding_eligibility(args.fasta, args.min_contig_len)
        kept = filter_fasta_by_length(args.fasta, filt_fa, args.min_contig_len)
        if kept != assessment["eligible_contigs"]:
            raise RuntimeError("Internal FASTA eligibility count mismatch")
        if kept < 2:
            unscaffolded = os.path.join(args.outdir, "unscaffolded_bin.fa")
            shutil.copy2(args.fasta, unscaffolded)
            write_scaffolding_status(
                os.path.join(args.outdir, "scaffolding_status.tsv"),
                args.fasta,
                "skipped",
                assessment,
            )
            logger.warning(
                "Skipping scaffolding: %d of %d contigs are at least %d bp; "
                "at least two eligible contigs are required",
                kept,
                assessment["total_contigs"],
                args.min_contig_len,
            )
            if not args.keep_temp:
                cleanup_intermediates(intermediate_dir)
            sys.exit(0)
        # Step 2: index
        if not os.path.exists(filt_fa + ".fai"):
            run_command(["samtools", "faidx", filt_fa])
        # Step 3: alignment or use BAM
        if args.bam:
            validate_bam_reference(args.bam, filt_fa)
            bam = os.path.join(intermediate_dir, "sorted_for_yahs.bam")
            run_command(["samtools", "sort", "-n", "-@", str(args.t), "-m", args.sort_memory,
                         args.bam, "-o", bam])
        else:
            aln_dir = os.path.join(intermediate_dir, "initial_alignment"); make_dir(aln_dir)
            run_command(alignment_command(args, filt_fa, aln_dir))
            bam = os.path.join(aln_dir, "sorted_map.bam")
        # Step 4: YaHS
        yahs_pref = os.path.join(intermediate_dir, "yahs", "scaffold")
        make_dir(os.path.dirname(yahs_pref))
        yahs_cmd = [
            "yahs",
            filt_fa,
            bam,
            "-o", yahs_pref,
            "-q", str(args.yahs_min_mapq),
            "-l", str(args.yahs_min_contig_len),
            "-R", str(args.yahs_rounds),
        ]
        if args.yahs_resolutions:
            yahs_cmd.extend(["-r", args.yahs_resolutions])
        if args.yahs_no_contig_ec:
            yahs_cmd.append("--no-contig-ec")
        if args.yahs_no_scaffold_ec:
            yahs_cmd.append("--no-scaffold-ec")
        if args.yahs_no_mem_check:
            yahs_cmd.append("--no-mem-check")
        if args.yahs_extra_args:
            yahs_cmd.extend(shlex.split(args.yahs_extra_args))
        run_command(yahs_cmd)
        yahs_genome = yahs_pref + "_scaffolds_final.fa"
        if not os.path.exists(yahs_genome): sys.exit(1)
        genome = os.path.join(args.outdir, "scaffolded_bin.fa")
        shutil.copy2(yahs_genome, genome)
        yahs_agp = yahs_pref + "_scaffolds_final.agp"
        if os.path.exists(yahs_agp):
            shutil.copy2(
                yahs_agp,
                os.path.join(args.outdir, "scaffolded_bin.agp"),
            )
        # Step 5: segment
        seg_len = args.r if args.r else 10000
        seg_fa = os.path.join(intermediate_dir, "segmented_scaffolds.fa")
        segmap = cut_scaffolds_to_segments(genome, seg_fa, seg_len)
        # Step 6: realign
        real_dir = os.path.join(intermediate_dir, "realign_to_scaffold"); make_dir(real_dir)
        run_command(alignment_command(args, seg_fa, real_dir))
        real_bam = os.path.join(real_dir, "sorted_map.bam")
        # Step 7: MetaCC
        metacc_dir = os.path.join(intermediate_dir, "metacc"); make_dir(metacc_dir)
        run_command([
            sys.executable,
            os.path.join(args.p, "5_contact", "raw_contact.py"),
            real_bam,
            ",".join(args.enzyme),
            seg_fa,
            metacc_dir,
            str(args.metacc_min_mapq),
            str(args.metacc_min_len),
            str(args.metacc_min_match),
            str(args.metacc_min_signal),
        ])
        # Step 8: normalize
        contig_csv = os.path.join(metacc_dir, "contig_info.csv")
        raw_contact = scisp.load_npz(
            os.path.join(metacc_dir, "Raw_contact_matrix.npz")
        ).tocsr()
        df = pd.read_csv(contig_csv)
        raw_contact, df = filter_metacc_inputs(raw_contact, df)
        scisp.save_npz(os.path.join(metacc_dir, "Raw_contact_matrix.npz"), raw_contact)
        df.to_csv(contig_csv, index=False)
        raw_contact.setdiag(0)
        raw_contact.eliminate_zeros()
        norm_res = normcc(contig_csv)
        hz = NormCCMap(metacc_dir, df, raw_contact, norm_res, thres=args.normcc_thres)
        scisp.save_npz(os.path.join(metacc_dir, "Normalized_contact_matrix.npz"), hz.seq_map.tocsr())
        with gzip.open(os.path.join(metacc_dir, "NormCC_normalized_contact.gz"), "wb") as f: pickle.dump(hz, f)
        # Step 9: mapping
        final_bins = os.path.join(intermediate_dir, "final_bins.p.gz")
        create_scaffold_mapping(segmap, final_bins)
        # Step 10: heatmap
        figs = os.path.join(args.outdir, "figures"); make_dir(figs)
        run_command([
            sys.executable,
            os.path.join(args.p, "8_scaffolding", "heatmap.py"),
            "--contact-map", os.path.join(metacc_dir, "Normalized_contact_matrix.npz"),
            "--contig-info", contig_csv,
            "--BIN", final_bins,
            "--OUTDIR", figs,
            "--max_image", str(args.heatmap_max_image),
        ])
        # Step 11: metrics
        metrics = calculating_metrics(
            filt_fa,
            genome,
            args,
            bam,
            real_bam,
            os.path.join(args.outdir, "quality"),
        )
        save_metrics(metrics, os.path.join(args.outdir, "scaffolding_metrics.txt"))
        assessment["reason"] = ""
        write_scaffolding_status(
            os.path.join(args.outdir, "scaffolding_status.tsv"),
            args.fasta,
            "completed",
            assessment,
        )
        if not args.keep_temp:
            cleanup_intermediates(intermediate_dir)
        logger.info("Scaffolding complete.")
    except Exception as e:
        logger.error(str(e)); sys.exit(1)
