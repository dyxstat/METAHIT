#!/usr/bin/env bash
# 10_MGE.sh - MGE detection and host interaction analysis
# Steps:
#   1. Run geNomad on combined contigs
#   1.5 Remove proviruses and make first report
#   2   Run CheckV QC on viral contigs, keep HQ/MQ only, make second report
#   3   Use normalized contact matrix to find host–MGE interactions (HQ vs MQ)
#   4   Generate final combined summary

set -euo pipefail

usage() {
    cat << 'EOF'
Usage: 10_MGE.sh [OPTIONS]

Required:
  -p, --metahit-path PATH      Path to metahit installation directory
  --combined PATH              Combined contigs FASTA (with bin/unmapped prefixes)
  --contact PATH               Normalized contact matrix (npz)
  --outdir PATH                Output directory

Optional:
  -t, --threads INT            Threads (default: 4)
EOF
}

THREADS=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--metahit-path) METAHIT_PATH="$2"; shift 2;;
    --combined)        COMBINED="$2"; shift 2;;
    --contact)         CONTACT_MATRIX="$2"; shift 2;;
    --outdir)          OUTDIR="$2"; shift 2;;
    -t|--threads)      THREADS="$2"; shift 2;;
    -h|--help)         usage; exit 0;;
    *) echo "Error: Unknown parameter: $1"; usage; exit 1;;
  esac
done

[[ -f "$COMBINED" ]] || { echo "Error: combined FASTA not found: $COMBINED"; exit 1; }
[[ -f "$CONTACT_MATRIX" ]] || { echo "Error: contact matrix not found: $CONTACT_MATRIX"; exit 1; }
[[ -n "${OUTDIR:-}" ]] || { echo "Error: --outdir is required"; exit 1; }

echo "[INFO] Output dir: $OUTDIR"
mkdir -p "$OUTDIR/genomad_output"

################################################################################
# STEP 1: Run geNomad
################################################################################
echo "[INFO] ===== STEP 1: Run geNomad ====="
eval "$(conda shell.bash hook)"
conda activate genomad

GENOMAD_DB="$METAHIT_PATH/databases/genomad_db"
[[ -f "$GENOMAD_DB/version.txt" ]] || { echo "[ERROR] geNomad DB missing at $GENOMAD_DB"; exit 1; }

genomad end-to-end \
  --cleanup \
  --splits 8 \
  --threads "$THREADS" \
  "$COMBINED" \
  "$OUTDIR/genomad_output" \
  "$GENOMAD_DB"

conda activate metahit_env
echo "[INFO] geNomad completed."

ASSEMBLY_BASE=$(basename "$COMBINED" .fa)
ASSEMBLY_BASE=$(basename "$ASSEMBLY_BASE" .fasta)
VIRUS_SUMMARY="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_virus_summary.tsv"
VIRAL_CONTIGS_RAW="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_virus.fna"
PLASMID_CONTIGS_RAW="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_plasmid.fna"

################################################################################
# STEP 1.5: Remove proviruses
################################################################################
echo "[INFO] ===== STEP 1.5: Remove proviruses ====="
FILTERED_VIRAL="$OUTDIR/checkv_output/virus_no_provirus.fna"
mkdir -p "$OUTDIR/checkv_output/virus"

set -a
VIRUS_SUMMARY="$VIRUS_SUMMARY"
VIRAL_CONTIGS_RAW="$VIRAL_CONTIGS_RAW"
FILTERED_VIRAL="$FILTERED_VIRAL"
set +a

python3 <<'EOF'
import pandas as pd
from Bio import SeqIO
import os

virus_summary = os.environ["VIRUS_SUMMARY"]
viral_fasta   = os.environ["VIRAL_CONTIGS_RAW"]
filtered_out  = os.environ["FILTERED_VIRAL"]

df = pd.read_csv(virus_summary, sep="\t")
keep_ids = set(df.loc[df["topology"] != "Provirus", "seq_name"])

with open(filtered_out, "w") as out_f:
    for rec in SeqIO.parse(viral_fasta, "fasta"):
        if rec.id in keep_ids:
            SeqIO.write(rec, out_f, "fasta")
EOF

################################################################################
# STEP 2: Run CheckV (only viruses)
################################################################################
echo "[INFO] ===== STEP 2: Run CheckV ====="
CHECKV_DB="$METAHIT_PATH/databases/checkv_db/checkv-db-v1.5"

conda activate checkv_env
checkv end_to_end "$FILTERED_VIRAL" "$OUTDIR/checkv_output/virus" -t "$THREADS" -d "$CHECKV_DB"
conda activate metahit_env

VIRAL_QC="$OUTDIR/checkv_output/virus_qc.fna"

################################################################################
# STEP 1.5 + STEP 2 + STEP 3 Reports in one summary
################################################################################
echo "[INFO] ===== Reporting ====="

set -a
OUTDIR="$OUTDIR"
VIRUS_SUMMARY="$VIRUS_SUMMARY"
VIRAL_CONTIGS_RAW="$VIRAL_CONTIGS_RAW"
FILTERED_VIRAL="$FILTERED_VIRAL"
PLASMID_CONTIGS_RAW="$PLASMID_CONTIGS_RAW"
VIRAL_QC="$VIRAL_QC"
CONTACT_MATRIX="$CONTACT_MATRIX"
COMBINED="$COMBINED"
set +a

python3 <<'EOF'
import os
import pandas as pd
from Bio import SeqIO
from scipy.sparse import load_npz

outdir        = os.environ["OUTDIR"]
virus_summary = os.environ["VIRUS_SUMMARY"]
filtered_viral= os.environ["FILTERED_VIRAL"]
plasmid_raw   = os.environ["PLASMID_CONTIGS_RAW"]
viral_qc      = os.environ["VIRAL_QC"]
contact_file  = os.environ["CONTACT_MATRIX"]
combined      = os.environ["COMBINED"]

summary_file = f"{outdir}/checkv_output/MGE_summary.txt"

def count_bin_unmapped(fasta):
    bins = sum(1 for rec in SeqIO.parse(fasta, "fasta") if rec.id.startswith("bin"))
    total = sum(1 for _ in SeqIO.parse(fasta, "fasta"))
    unmapped = total - bins
    return total, bins, unmapped

with open(summary_file,"w") as fh:
    fh.write("=== STEP 1.5 Report (After Provirus Removal) ===\n")
    df = pd.read_csv(virus_summary, sep="\t")
    provirus_count = (df["topology"] == "Provirus").sum()
    virus_total, virus_bins, virus_unmapped = count_bin_unmapped(filtered_viral)
    plasmid_total, plasmid_bins, plasmid_unmapped = count_bin_unmapped(plasmid_raw)
    fh.write(f"Virus contigs (after provirus removal): {virus_total}\n")
    fh.write(f"  from bins: {virus_bins}\n")
    fh.write(f"  from unmapped: {virus_unmapped}\n")
    fh.write(f"Provirus removed: {provirus_count}\n")
    fh.write(f"Plasmid contigs (raw): {plasmid_total}\n")
    fh.write(f"  from bins: {plasmid_bins}\n")
    fh.write(f"  from unmapped: {plasmid_unmapped}\n\n")

    fh.write("=== STEP 2 Report (CheckV Quality, Viruses Only) ===\n")
    virus_df = pd.read_csv(f"{outdir}/checkv_output/virus/quality_summary.tsv", sep="\t")
    virus_keep = virus_df[virus_df["checkv_quality"].isin(["High-quality","Medium-quality"])]
    hq_count = (virus_keep["checkv_quality"] == "High-quality").sum()
    mq_count = (virus_keep["checkv_quality"] == "Medium-quality").sum()
    fh.write(f"Virus HQ: {hq_count}\n")
    fh.write(f"Virus MQ: {mq_count}\n")
    fh.write(f"Virus HQ+MQ total: {len(virus_keep)}\n\n")

    # Export QC FASTA for HQ+MQ viruses
    ids_virus = set(virus_keep["contig_id"])
    with open(viral_qc,"w") as out_f:
        for rec in SeqIO.parse(filtered_viral,"fasta"):
            if rec.id in ids_virus:
                SeqIO.write(rec,out_f,"fasta")

    fh.write("=== STEP 3 Report (Host–MGE Contacts) ===\n")
    mat = load_npz(contact_file).tocoo()
    contigs = [line[1:].split()[0] for line in open(combined) if line.startswith(">")]
    viral_contigs = set([line[1:].split()[0] for line in open(viral_qc) if line.startswith(">")])
    all_mge = viral_contigs
    host_contigs = {c for c in contigs if (c.startswith("bin") and c not in all_mge)}
    quality_map = dict(zip(virus_df["contig_id"], virus_df["checkv_quality"]))
    hq_set = {c for c in viral_contigs if quality_map.get(c)=="High-quality"}
    mq_set = {c for c in viral_contigs if quality_map.get(c)=="Medium-quality"}
    def count_contacts(mge_set):
        count=0
        for r,c,v in zip(mat.row, mat.col, mat.data):
            if v<=0: continue
            a,b=contigs[r], contigs[c]
            if (a in mge_set and b in host_contigs) or (b in mge_set and a in host_contigs):
                count+=1
        return count
    fh.write(f"HQ viral contigs: {len(hq_set)}\n")
    fh.write(f"MQ viral contigs: {len(mq_set)}\n")
    fh.write(f"HQ–host contacts: {count_contacts(hq_set)}\n")
    fh.write(f"MQ–host contacts: {count_contacts(mq_set)}\n")
EOF

echo "[INFO] ===== MGE analysis finished, summary in checkv_output/MGE_summary.txt ====="
