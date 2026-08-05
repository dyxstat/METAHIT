#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  7_reassembly.sh -p <METAHICT_PATH> --bin <BIN> --assembly <ASSEMBLY> --hic1 <HIC1> --hic2 <HIC2> --sg1 <SG1> --sg2 <SG2> --bam <BAM> --outdir <OUTDIR> [options]

Required arguments:
  -p PATH                    METAHICT modules path.
  --bin PATH                 Input bin directory.
  --assembly PATH            Assembly FASTA used to create the contig sort list.
  --hic1 PATH                Hi-C forward FASTQ.gz.
  --hic2 PATH                Hi-C reverse FASTQ.gz.
  --sg1 PATH                 Shotgun forward FASTQ.gz.
  --sg2 PATH                 Shotgun reverse FASTQ.gz.
  --bam PATH                 Hi-C BAM mapped to the assembly.
  --outdir PATH              Output directory.

Options:
  -t, --threads INT          Number of CPU threads (default=80).
  -m, --memory INT           Memory in GB (default=80% of available system memory).
  --cutoff-quantile FLOAT    Quantile of the short-insert component used for selection (default=0.95).
  --top-k INT                Number of longest contigs used for EM fitting (default=100).
  --min-mapq INT             Minimum MAPQ for pair-level insert-size extraction (default=30).
  --min-match-len INT        Minimum aligned match length for pair-level extraction (default=30).
  --exclude-duplicates       Exclude duplicate-marked alignments.
  --write-nonselected-hic    Also write non-selected Hi-C FASTQ files.
  --min-contig-len INT       Minimum reassembled contig length retained (default=500).
  --strict-cut-off INT       Strict read-recruitment SNP cutoff (default=2).
  --permissive-cut-off INT   Permissive read-recruitment SNP cutoff (default=5).
  --contamination-penalty FLOAT
                              Penalty used in completeness - penalty * contamination (default=5).
  --skip-checkm2             Skip CheckM2 and final best-bin selection.
  --checkm2_db PATH          CheckM2 database path.
  --tmp-dir PATH             Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp).
  --spades-mode TEXT         SPAdes mode: careful or none (default=careful).
  --spades-phred-offset INT  SPAdes phred offset (default=empty).
  --spades-extra-args TEXT   Additional options passed directly to SPAdes (default=empty).
  --skip-residual-assembly   Do not assemble residual unmapped reads.
  --keep-temp                Keep successful SPAdes and CheckM2 temporary directories.
  --print-defaults           Print default values and exit.
  -h, --help                 Show this help and exit.
EOF
}

print_defaults() {
    cat <<'EOF'
Module 7 reassembly defaults:
threads=80
memory=80% of available system memory
cutoff_quantile=0.95
top_k=100
min_mapq=30
min_match_len=30
exclude_duplicates=false
write_nonselected_hic=false
min_contig_len=500
strict_cut_off=2
permissive_cut_off=5
contamination_penalty=5
skip_checkm2=false
checkm2_db=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
spades_mode=careful
spades_phred_offset=
spades_extra_args=
skip_residual_assembly=false
keep_temp=false
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

source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"
METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "$METAHICT_RESOURCE_OUTDIR" ]]; then
    metahict_resource_start "$METAHICT_RESOURCE_OUTDIR" "$(basename "$0")"
fi

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

CHECKM2_DB=""
THREADS=80
MEMORY=""
CUTOFF_QUANTILE="0.95"
TOP_K="100"
MIN_MAPQ="30"
MIN_MATCH_LEN="30"
EXCLUDE_DUPLICATES=false
WRITE_NONSELECTED_HIC=false
MIN_CONTIG_LEN="500"
STRICT_CUT_OFF="2"
PERMISSIVE_CUT_OFF="5"
CONTAMINATION_PENALTY="5"
SKIP_CHECKM2=false
TMP_DIR=""
SPADES_MODE="careful"
SPADES_PHRED_OFFSET=""
SPADES_EXTRA_ARGS=""
SKIP_RESIDUAL_ASSEMBLY=false
KEEP_TEMP=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p) PATH_DIR="$2"; shift 2 ;;
        --bin) BIN="$2"; shift 2 ;;
        --assembly) ASSEMBLY="$2"; shift 2 ;;
        --hic1) HIC1="$2"; shift 2 ;;
        --hic2) HIC2="$2"; shift 2 ;;
        --sg1) SG1="$2"; shift 2 ;;
        --sg2) SG2="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -m|--memory) MEMORY="$2"; shift 2 ;;
        --cutoff-quantile) CUTOFF_QUANTILE="$2"; shift 2 ;;
        --top-k|--top_k) TOP_K="$2"; shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-match-len) MIN_MATCH_LEN="$2"; shift 2 ;;
        --exclude-duplicates) EXCLUDE_DUPLICATES=true; shift 1 ;;
        --write-nonselected-hic) WRITE_NONSELECTED_HIC=true; shift 1 ;;
        --min-contig-len) MIN_CONTIG_LEN="$2"; shift 2 ;;
        --strict-cut-off) STRICT_CUT_OFF="$2"; shift 2 ;;
        --permissive-cut-off) PERMISSIVE_CUT_OFF="$2"; shift 2 ;;
        --contamination-penalty) CONTAMINATION_PENALTY="$2"; shift 2 ;;
        --skip-checkm2) SKIP_CHECKM2=true; shift 1 ;;
        --checkm2_db) CHECKM2_DB="$2"; shift 2 ;;
        --tmp-dir) TMP_DIR="$2"; shift 2 ;;
        --spades-mode) SPADES_MODE="$2"; shift 2 ;;
        --spades-phred-offset) SPADES_PHRED_OFFSET="$2"; shift 2 ;;
        --spades-extra-args) SPADES_EXTRA_ARGS="$2"; shift 2 ;;
        --skip-residual-assembly) SKIP_RESIDUAL_ASSEMBLY=true; shift 1 ;;
        --keep-temp) KEEP_TEMP=true; shift 1 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "${PATH_DIR:-}" || -z "${BIN:-}" || -z "${ASSEMBLY:-}" || -z "${HIC1:-}" || -z "${HIC2:-}" || -z "${SG1:-}" || -z "${SG2:-}" || -z "${BAM:-}" || -z "${OUTDIR:-}" ]]; then
    echo "Error: Missing required arguments"
    usage
    exit 1
fi

if [[ -z "$MEMORY" ]]; then
    MEMORY=$(awk '/MemTotal/ {printf "%d", ($2 / 1024 / 1024) * 0.8}' /proc/meminfo)
fi

mkdir -p "$OUTDIR"

ASSEMBLY_WORKDIR="${OUTDIR}/input_assembly"
mkdir -p "$ASSEMBLY_WORKDIR"
WORK_ASSEMBLY="${ASSEMBLY_WORKDIR}/$(basename "$ASSEMBLY")"
cp "$ASSEMBLY" "$WORK_ASSEMBLY"

CONTIG_SORT_LIST="${OUTDIR}/contig_sort.txt"
echo "Creating contig sort list: $CONTIG_SORT_LIST"
if [[ ! -f "${WORK_ASSEMBLY}.fai" ]]; then
    samtools faidx "$WORK_ASSEMBLY"
fi
sort -k2,2nr "${WORK_ASSEMBLY}.fai" | cut -f1 > "$CONTIG_SORT_LIST"

REASSEMBLY_SCRIPT="${PATH_DIR}/7_reassembly/scripts/7_reassembly.py"

CMD=(python "$REASSEMBLY_SCRIPT"
    --bin "$BIN"
    --hic1 "$HIC1"
    --hic2 "$HIC2"
    --sg1 "$SG1"
    --sg2 "$SG2"
    --bam "$BAM"
    --outdir "$OUTDIR"
    -p "$PATH_DIR"
    -c "$CONTIG_SORT_LIST"
    -t "$THREADS"
    -m "$MEMORY"
    --top_k "$TOP_K"
    --min-mapq "$MIN_MAPQ"
    --min-match-len "$MIN_MATCH_LEN"
    --cutoff-quantile "$CUTOFF_QUANTILE"
    --min-contig-len "$MIN_CONTIG_LEN"
    --strict-cut-off "$STRICT_CUT_OFF"
    --permissive-cut-off "$PERMISSIVE_CUT_OFF"
    --contamination-penalty "$CONTAMINATION_PENALTY"
    --spades-mode "$SPADES_MODE"
    # Current 7_reassembly.py uses action=store_false for this flag, so passing
    # it forces creation of a run-local name-sorted BAM.
    --bam-name-sorted
)

if [[ "$EXCLUDE_DUPLICATES" = true ]]; then
    CMD+=(--exclude-duplicates)
fi

if [[ "$WRITE_NONSELECTED_HIC" = true ]]; then
    CMD+=(--write-nonselected-hic)
fi

if [[ "$SKIP_CHECKM2" = true ]]; then
    CMD+=(--skip-checkm2)
fi

if [[ -n "$CHECKM2_DB" ]]; then
    CMD+=(--checkm2_db "$CHECKM2_DB")
fi

if [[ -n "$TMP_DIR" ]]; then
    CMD+=(--tmp-dir "$TMP_DIR")
fi

if [[ -n "$SPADES_PHRED_OFFSET" ]]; then
    CMD+=(--spades-phred-offset "$SPADES_PHRED_OFFSET")
fi

if [[ -n "$SPADES_EXTRA_ARGS" ]]; then
    CMD+=("--spades-extra-args=${SPADES_EXTRA_ARGS}")
fi

if [[ "$SKIP_RESIDUAL_ASSEMBLY" = true ]]; then
    CMD+=(--skip-residual-assembly)
fi

if [[ "$KEEP_TEMP" = true ]]; then
    CMD+=(--keep-temp)
fi

echo "Running reassembly with command:"
printf '%q ' "${CMD[@]}"
echo
"${CMD[@]}"

echo "Reassembly completed successfully."
