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
import warnings
import gzip
import pickle
import pandas as pd
import scipy.sparse as scisp
from pandas.errors import SettingWithCopyWarning
warnings.filterwarnings("ignore", category=SettingWithCopyWarning)

from raw_contact_both import ContactMatrix
from MetaCC.Script.normalized_contact import NormCCMap
from MetaCC.Script.utils import save_object, make_dir
from MetaCC.Script.normcc import normcc

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
    cmd = f"conda run -n checkm2 checkm2 predict --threads {threads} --input {fasta} --output-directory {outdir}"
    try:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)
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


def ensure_bam_indexed(bam):
    idx = bam + ".bai"
    if not os.path.exists(idx):
        logger.info(f"Indexing BAM {bam}")
        if os.system(f"samtools index {bam}") != 0:
            logger.error(f"Failed to index {bam}")
            return False
    return True


def create_coordinate_sorted_bam(bam_in, bam_out, threads=1):
    if os.path.exists(bam_out):
        logger.info(f"Sorted BAM exists: {bam_out}")
        return True
    logger.info(f"Sorting BAM {bam_in}")
    if os.system(f"samtools sort -@ {threads} {bam_in} -o {bam_out}") != 0:
        logger.error(f"Failed to sort {bam_in}")
        return False
    return True


def calculate_hic_enrichment_contigs(bam, fasta):
    import pysam
    if not ensure_bam_indexed(bam):
        return None
    contigs = {l[1:].strip().split()[0] for l in open(fasta) if l.startswith(">")}
    within, across = 0, 0
    with pysam.AlignmentFile(bam, "rb") as bf:
        for r in bf.fetch():
            if r.is_paired and not r.is_unmapped and not r.mate_is_unmapped:
                c1, c2 = r.reference_name, r.next_reference_name
                if c1 in contigs and c2 in contigs:
                    if c1 == c2: within += 1
                    else: across += 1
    return within / across if across else (float("inf") if within else 0)


def calculate_hic_enrichment_scaffolds(bam, fasta):
    import pysam
    if not ensure_bam_indexed(bam):
        return None
    scaffs = {l[1:].strip().split()[0] for l in open(fasta) if l.startswith(">")}
    def sname(seg): return seg.split("_seg_")[0] if "_seg_" in seg else seg
    within, across = 0, 0
    with pysam.AlignmentFile(bam, "rb") as bf:
        for r in bf.fetch():
            if r.is_paired and not r.is_unmapped and not r.mate_is_unmapped:
                s1, s2 = sname(r.reference_name), sname(r.next_reference_name)
                if s1 in scaffs and s2 in scaffs:
                    if s1 == s2: within += 1
                    else: across += 1
    return within / across if across else (float("inf") if within else 0)


def calculating_metrics(filt_fa, genome, cm, segmap, args, outdir):
    logger.info("Calculating quality metrics")
    n50_o, l50_o = calculate_n50_l50(filt_fa)
    n50_s, l50_s = calculate_n50_l50(genome)
    c_o = sum(1 for l in open(filt_fa) if l.startswith(">"))
    c_s = sum(1 for l in open(genome) if l.startswith(">"))
    d1, d2 = os.path.join(outdir, "checkm2_original"), os.path.join(outdir, "checkm2_scaffolded")
    make_dir(d1); make_dir(d2)
    t = int(args.t) if args.t else 1
    db = getattr(args, "checkm2_db", None)
    comp_o, cont_o = run_checkm2(filt_fa, d1, t, db)
    comp_s, cont_s = run_checkm2(genome, d2, t, db)
    enr_o, enr_s = None, None
    bam_o = os.path.join(outdir, "sorted_for_yahs.bam")
    bam_s = os.path.join(outdir, "realign_to_scaffold", "sorted_map.bam")
    if os.path.exists(bam_o): enr_o = calculate_hic_enrichment_contigs(bam_o, filt_fa)
    if os.path.exists(bam_s): enr_s = calculate_hic_enrichment_scaffolds(bam_s, genome)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="Base metahit path")
    parser.add_argument("--fasta", help="Bin FASTA")
    parser.add_argument("--bam", help="Optional BAM")
    parser.add_argument("--enzyme", help="Restriction enzyme list")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--hic1", help="Hi-C forward read")
    parser.add_argument("--hic2", help="Hi-C reverse read")
    parser.add_argument("-t", help="Threads")
    parser.add_argument("-m", help="Memory")
    parser.add_argument("-r", help="Resolution (default 10kb)")
    parser.add_argument("--checkm2_db", help="Path to custom CheckM2 database")
    args = parser.parse_args()

    make_dir(args.outdir)
    log_file = os.path.join(args.outdir, "scaffolding.log")
    logging.basicConfig(level=logging.INFO, format="%(levelname)-8s | %(asctime)s | %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="a")])
    global logger
    logger = logging.getLogger("main")

    try:
        # Step 1: filter
        filt_fa = os.path.join(args.outdir, "filtered_bin.fa")
        kept = filter_fasta_by_length(args.fasta, filt_fa, 5000)
        if kept == 0:
            logger.error("No contigs >5kb"); sys.exit(1)
        # Step 2: index
        if not os.path.exists(filt_fa + ".fai"):
            os.system(f"samtools faidx {filt_fa}")
        # Step 3: alignment or use BAM
        if args.bam:
            bam = os.path.join(args.outdir, "sorted_for_yahs.bam")
            os.system(f"samtools sort -n {args.bam} -o {bam}")
        else:
            aln_dir = os.path.join(args.outdir, "initial_alignment"); make_dir(aln_dir)
            aln_cmd = (f"bash {os.path.dirname(__file__)}/alignment.sh -p {args.p} "
                       f"-r {filt_fa} -1 {args.hic1} -2 {args.hic2} -o {aln_dir} -t {args.t}")
            if os.system(aln_cmd) != 0: sys.exit(1)
            bam = os.path.join(aln_dir, "sorted_map.bam")
        # Step 4: YaHS
        yahs_pref = os.path.join(args.outdir, "yahs", "scaffold")
        make_dir(os.path.dirname(yahs_pref))
        if os.system(f"yahs {filt_fa} {bam} -o {yahs_pref}") != 0: sys.exit(1)
        genome = yahs_pref + "_scaffolds_final.fa"
        if not os.path.exists(genome): sys.exit(1)
        # Step 5: segment
        seg_len = int(args.r) if args.r else 10000
        seg_fa = os.path.join(args.outdir, "segmented_scaffolds.fa")
        segmap = cut_scaffolds_to_segments(genome, seg_fa, seg_len)
        # Step 6: realign
        real_dir = os.path.join(args.outdir, "realign_to_scaffold"); make_dir(real_dir)
        real_cmd = (f"bash {os.path.dirname(__file__)}/alignment.sh -p {args.p} "
                    f"-r {seg_fa} -1 {args.hic1} -2 {args.hic2} -o {real_dir} -t {args.t}")
        if os.system(real_cmd) != 0: sys.exit(1)
        real_bam = os.path.join(real_dir, "sorted_map.bam")
        # Step 7: MetaCC
        metacc_dir = os.path.join(args.outdir, "metacc"); make_dir(metacc_dir)
        cm = ContactMatrix(real_bam, args.enzyme.split(",") if args.enzyme else ["Sau3AI", "MluCI"],
                           seg_fa, args.outdir, metacc_dir, None)
        cm_path = os.path.join(metacc_dir, "contact_map.p.gz")
        with gzip.open(cm_path, "wb") as f: pickle.dump(cm, f)
        # Step 8: normalize
        contig_csv = os.path.join(metacc_dir, "tmp", "contig_info_metacc.csv")
        norm_res = normcc(contig_csv)
        df = pd.read_csv(contig_csv)
        hz = NormCCMap(metacc_dir, df, cm.seq_map_metacc, norm_res, thres=0.05)
        scisp.save_npz(os.path.join(metacc_dir, "Normalized_contact_matrix.npz"), hz.seq_map.tocsr())
        with gzip.open(os.path.join(metacc_dir, "NormCC_normalized_contact.gz"), "wb") as f: pickle.dump(hz, f)
        # Step 9: mapping
        final_bins = os.path.join(args.outdir, "final_bins.p.gz")
        create_scaffold_mapping(segmap, final_bins)
        # Step 10: heatmap
        figs = os.path.join(args.outdir, "figures"); make_dir(figs)
        heat_cmd = (f"python {args.p}/8_scaffolding/scripts/heatmap.py "
                    f"--contact-map {os.path.join(metacc_dir, 'Normalized_contact_matrix.npz')} "
                    f"--ORDER {cm_path} --BIN {final_bins} --OUTDIR {figs}")
        if os.system(heat_cmd) != 0: sys.exit(1)
        # Step 11: metrics
        metrics = calculating_metrics(filt_fa, genome, cm, segmap, args, args.outdir)
        save_metrics(metrics, os.path.join(args.outdir, "scaffolding_metrics.txt"))
        logger.info("Scaffolding complete.")
    except Exception as e:
        logger.error(str(e)); sys.exit(1)
