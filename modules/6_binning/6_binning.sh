#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage:
  6_binning.sh <FASTA> <BAM> <OUTDIR> <PROJECT_PATH> [options]

Required arguments:
  <FASTA>                 Assembled contigs FASTA.
  <BAM>                   Hi-C read alignments to the contigs.
  <OUTDIR>                Output directory for binning results.
  <PROJECT_PATH>          METAHICT project path.

Options:
  -t, --threads INT       Number of CPU threads (default=80).
  --checkm_db PATH        CheckM database path.
  --enzyme STR            Comma-separated restriction enzymes (default=Sau3AI,MluCI).
  --metacc-min-len INT    Minimum contig length for MetaCC contact generation (default=1000).
  --metacc-min-signal INT Minimum Hi-C signal for MetaCC contact generation (default=2).
  --metacc-min-mapq INT   Minimum MAPQ for MetaCC contact generation (default=30).
  --metacc-min-match INT  Minimum aligned match length for MetaCC contact generation (default=30).
  --metacc-min-binsize INT Minimum MetaCC output bin size (default=150000).
  --normcc-thres FLOAT    Fraction of NormCC-normalized contacts discarded (default=0.05).
  --bin3c-min-len INT     Minimum contig length for bin3C-compatible contact generation (default=1000).
  --bin3c-min-signal INT  Minimum Hi-C signal for bin3C-compatible contact generation (default=5).
  --bin3c-min-mapq INT    Minimum MAPQ for bin3C-compatible contact generation (default=60).
  --bin3c-min-match INT   Minimum aligned match length for bin3C-compatible contact generation (default=10).
  --bin3c-min-extent INT  Minimum bin3C cluster extent (default=50000).
  --min-completeness FLOAT Final METAHICT minimum completeness (default=50).
  --max-contamination FLOAT Final METAHICT maximum contamination (default=10).
  --contamination-penalty FLOAT Penalty used in completeness - penalty * contamination (default=5).
  --min-input-bin-size INT Minimum input bin FASTA file size before refinement (default=50000 bytes).
  --max-input-bin-size INT Maximum input bin FASTA file size before refinement (default=20000000 bytes).
  --binning-refiner-min-size INT Minimum refined bin size for Binning_refiner (default=524288 bp).
  --tmp-dir PATH          Temporary directory root for CheckM2 working files (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp).
  --keep-temp             Keep successful CheckM2 temporary directories.
  --no-fasta              Do not write bin3C cluster FASTA files.
  --no-report             Do not write bin3C cluster report.
  --no-spades             Input assembly was not produced by SPAdes.
  --only-large            Only write bin3C FASTA clusters larger than --bin3c-min-extent.
  --skip-checkm2          Skip CheckM2 during final bin refinement.
  --skip-refinement       Skip Binning_refiner combinations.
  --skip-consolidation    Skip final consolidation across bin sets.
  --keep-ambiguous        Keep ambiguous contigs in all bins.
  --remove-ambiguous      Remove ambiguous contigs from all bins.
  --seed INT              Random seed.
  --print-defaults        Print default values and exit.
  -h, --help              Show this help and exit.
EOF
}

print_defaults() {
    cat <<'EOF'
Module 6 binning defaults:
threads=80
enzyme=Sau3AI,MluCI
metacc_min_len=1000
metacc_min_signal=2
metacc_min_mapq=30
metacc_min_match=30
metacc_min_binsize=150000
normcc_thres=0.05
bin3c_min_len=1000
bin3c_min_signal=5
bin3c_min_mapq=60
bin3c_min_match=10
bin3c_min_extent=50000
min_completeness=50
max_contamination=10
contamination_penalty=5
min_input_bin_size=50000
max_input_bin_size=20000000
binning_refiner_min_size=524288
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
no_fasta=false
no_report=false
no_spades=false
only_large=false
keep_temp=false
skip_checkm2=false
skip_refinement=false
skip_consolidation=false
keep_ambiguous=false
remove_ambiguous=false
seed=
EOF
}

for arg in "$@"; do
    case "$arg" in
        -h|--help)
            usage
            exit 0
            ;;
        --print-defaults)
            print_defaults
            exit 0
            ;;
    esac
done

# Check for required arguments
if [ "$#" -lt 4 ]; then
    usage
    exit 1
fi

OUTDIR_FOR_CLEANUP="$3"
if [[ -d "$OUTDIR_FOR_CLEANUP" ]]; then
    echo "[INFO] Removing previous binning output: $OUTDIR_FOR_CLEANUP"
    rm -rf "$OUTDIR_FOR_CLEANUP"
fi

source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"
METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "$METAHICT_RESOURCE_OUTDIR" ]]; then
    metahict_resource_start "$METAHICT_RESOURCE_OUTDIR" "$(basename "$0")"
fi

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

FASTA=$1
BAM=$2
OUTDIR=$3
PROJECT_PATH=$4
shift 4

THREADS=80
CHECKM_DB_PATH=""
BINNING_ARGS=()
INTEGRATION_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --checkm_db)
            CHECKM_DB_PATH="$2"
            BINNING_ARGS+=("--checkm_db" "$2")
            shift 2
            ;;
        --enzyme|--metacc-min-len|--metacc-min-signal|--metacc-min-mapq|--metacc-min-match|--metacc-min-binsize|--bin3c-min-len|--bin3c-min-signal|--bin3c-min-mapq|--bin3c-min-match|--bin3c-min-extent|--seed)
            BINNING_ARGS+=("$1" "$2")
            shift 2
            ;;
        --normcc-thres)
            BINNING_ARGS+=("--thres" "$2")
            shift 2
            ;;
        --min-completeness|--max-contamination|--contamination-penalty|--min-input-bin-size|--max-input-bin-size|--binning-refiner-min-size|--tmp-dir)
            INTEGRATION_ARGS+=("$1" "$2")
            shift 2
            ;;
        --no-fasta|--no-report|--no-spades|--only-large)
            BINNING_ARGS+=("$1")
            shift
            ;;
        --skip-checkm2|--skip-refinement|--skip-consolidation|--keep-ambiguous|--remove-ambiguous|--keep-temp)
            INTEGRATION_ARGS+=("$1")
            shift
            ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

BINNING_ARGS+=("--threads" "$THREADS")
INTEGRATION_ARGS+=("--threads" "$THREADS")

# ---- Part 1: Binning ----
BINNING_SCRIPT="${PROJECT_PATH}/6_binning/scripts/6a_binning.py"

echo "[INFO] Running binning tools (MetaCC, bin3C, ImputeCC)..."
python "$BINNING_SCRIPT" --FASTA "$FASTA" --BAM "$BAM" --OUTDIR "$OUTDIR" "${BINNING_ARGS[@]}" || {
    echo "[ERROR] Binning failed."
    exit 1
}
echo "[INFO] Binning completed successfully."

# ---- Part 2: Run METAHICT bin integration ----
METAHICT_SCRIPT="${PROJECT_PATH}/6_binning/scripts/6b_integration.py"
OUTDIR_METACC="${OUTDIR}/metacc"
OUTDIR_BIN3C="${OUTDIR}/bin3c"
OUTDIR_IMPUTECC="${OUTDIR}/imputecc"
OUTDIR_METAHICT="${OUTDIR}/metahict"

mkdir -p "$OUTDIR_METAHICT"

echo "[INFO] Running METAHICT bin integration..."
python "$METAHICT_SCRIPT" "$OUTDIR_METACC" "$OUTDIR_BIN3C" "$OUTDIR_IMPUTECC" "$OUTDIR_METAHICT" "${INTEGRATION_ARGS[@]}" || {
    echo "[ERROR] METAHICT bin integration failed."
    exit 1
}
echo "[INFO] METAHICT bin integration completed successfully."

# ---- Part 3: Plot Heatmap of Final Bins ----
HEATMAP_SCRIPT="${PROJECT_PATH}/6_binning/scripts/heatmap.py"
CONTACT_MATRIX="${OUTDIR_METACC}/Normalized_contact_matrix.npz"
CONTACT_MAP="${OUTDIR_METACC}/contact_map.p.gz"
CLUSTERING="${OUTDIR_METAHICT}/final_bins.p.gz"
PLOT_OUTDIR="${OUTDIR_METAHICT}/figures"
mkdir -p "$PLOT_OUTDIR"

echo "[INFO] Plotting heatmap of METAHICT final bins..."
python "$HEATMAP_SCRIPT" \
    --contact-map "$CONTACT_MATRIX" \
    --ORDER "$CONTACT_MAP" \
    --BIN "$CLUSTERING" \
    --OUTDIR "$PLOT_OUTDIR" || {
        echo "[ERROR] Heatmap plotting failed."
        exit 1
}
echo "[INFO] Heatmap plotting completed successfully."
