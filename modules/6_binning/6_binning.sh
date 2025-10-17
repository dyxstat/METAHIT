#!/usr/bin/env bash
set -euo pipefail

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

# Check for required arguments
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <FASTA> <BAM> <OUTDIR_BASE> <PROJECT_PATH> [optional args]"
    echo "Example:"
    echo "  $0 final.contigs.fa ALL_MAP_SORTED.bam bin_refinement metahit [--threads 80]"
    exit 1
fi

FASTA=$1
BAM=$2
OUTDIR=$3
PROJECT_PATH=$4
shift 4

# ---- Part 1: Binning ----
BINNING_SCRIPT="${PROJECT_PATH}/6_binning/scripts/6a_binning.py"

echo "[INFO] Running binning tools (MetaCC, bin3C, ImputeCC)..."
python "$BINNING_SCRIPT" --FASTA "$FASTA" --BAM "$BAM" --OUTDIR "$OUTDIR" "$@" || {
    echo "[ERROR] Binning failed."
    exit 1
}
echo "[INFO] Binning completed successfully."

# ---- Part 2: Run METAHIT bin integration ----
METAHIT_SCRIPT="${PROJECT_PATH}/6_binning/scripts/6b_integration.py"
OUTDIR_METACC="${OUTDIR}/metacc"
OUTDIR_BIN3C="${OUTDIR}/bin3c"
OUTDIR_IMPUTECC="${OUTDIR}/imputecc"
OUTDIR_METAHIT="${OUTDIR}/metahit"

mkdir -p "$OUTDIR_METAHIT"

echo "[INFO] Running MetaHIT bin integration..."
python "$METAHIT_SCRIPT" "$OUTDIR_METACC" "$OUTDIR_BIN3C" "$OUTDIR_IMPUTECC" "$OUTDIR_METAHIT" "$@" || {
    echo "[ERROR] MetaHIT bin integration failed."
    exit 1
}
echo "[INFO] MetaHIT bin integration completed successfully."

# ---- Part 3: Plot Heatmap of Final Bins ----
HEATMAP_SCRIPT="${PROJECT_PATH}/6_binning/scripts/heatmap.py"
CONTACT_MATRIX="${OUTDIR_METACC}/Normalized_contact_matrix.npz"
CONTACT_MAP="${OUTDIR_METACC}/contact_map.p.gz"
CLUSTERING="${OUTDIR_METAHIT}/final_bins.p.gz"
PLOT_OUTDIR="${OUTDIR_METAHIT}/figures"
mkdir -p "$PLOT_OUTDIR"

echo "[INFO] Plotting heatmap of MetaHIT final bins..."
python "$HEATMAP_SCRIPT" \
    --contact-map "$CONTACT_MATRIX" \
    --ORDER "$CONTACT_MAP" \
    --BIN "$CLUSTERING" \
    --OUTDIR "$PLOT_OUTDIR" || {
        echo "[ERROR] Heatmap plotting failed."
        exit 1
}
echo "[INFO] Heatmap plotting completed successfully."
