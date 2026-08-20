#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"
METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "$METAHICT_RESOURCE_OUTDIR" ]]; then
    metahict_resource_start "$METAHICT_RESOURCE_OUTDIR" "$(basename "$0")"
fi

set -euo pipefail

default_tmp_root() {
    local candidate="${METAHICT_TMP_ROOT:-${TMPDIR:-/tmp}}"
    if [[ -d "$candidate" && -w "$candidate" ]]; then
        printf "%s\n" "$candidate"
        return
    fi
    if [[ -d /tmp && -w /tmp ]]; then
        printf "%s\n" "/tmp"
        return
    fi
    return 1
}

print_defaults() {
    cat <<'EOF'
Module 3 alignment defaults:
threads=80
bwa_options=-5SP
samtools_filter=-F 0x900
mapq=30
min_intra_dist=10000
min_match_len=30
sort_memory=1G
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
keep_sam=false
skip_metrics=false
EOF
}

usage() {
    cat <<'EOF'

Usage: metahict alignment [options] -r reference.fa -1 reads_1.fastq -2 reads_2.fastq -o output_dir
Options:
  -p STR               MetaHit project path or modules path
  -r, --reference STR  Reference FASTA (default=output/assembly/final_assembly.fasta)
  -1, --reads1 STR     First Hi-C FASTQ (default=output/readqc/hic/final_reads_1.fastq)
  -2, --reads2 STR     Second Hi-C FASTQ (default=output/readqc/hic/final_reads_2.fastq)
  -o, --output STR     Output directory (default=output/alignment)
  -t, --threads INT    Number of CPU threads (default=80)
  --bwa-options STR    BWA MEM options (default="-5SP")
  --samtools-filter STR
                       samtools view filter options for BAM output (default="-F 0x900")
  --mapq INT           Minimum MAPQ retained in BAM and used for metrics (default=30)
  --min-intra-dist INT Minimum same-contig distance for informative pairs (default=10000)
  --min-match-len INT  Minimum aligned nucleotide match length retained in BAM and used for metrics (default=30)
  --sort-memory STR    Memory per samtools sort thread (default=1G)
  --tmp-dir STR        Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)
  --keep-sam           Keep intermediate map.sam (default=false)
  --skip-metrics       Skip 3D ratio and informative-pair metric calculation (default=false)
  --print-defaults     Print default parameter values and exit
  -h, --help           Show this help message

EOF
}

calculate_hic_ratios() {
    local bam_path="$1"
    local output_dir="$2"
    local metrics_script="${RUN_TMP}/calculate_hic_ratios.py"

    echo "[INFO] Calculating 3D ratio and informative pairs ratio..."
    echo "[INFO] Metric filters: MAPQ >= ${MAPQ}; min intra-contig distance = ${MIN_INTRA_DIST}; min match length = ${MIN_MATCH_LEN}"

    cat > "$metrics_script" <<'PY'
import re
import sys

output_dir = sys.argv[1]
min_intra_dist = int(sys.argv[2])
min_match_len = int(sys.argv[3])

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")

def match_length(cigar):
    if cigar == "*":
        return 0
    return sum(int(n) for n, op in cigar_re.findall(cigar) if op in {"M", "=", "X"})

def passes_match_length(fields):
    return min_match_len <= 0 or match_length(fields[5]) >= min_match_len

def is_informative_pair(a_fields, b_fields, check_distance=True):
    flag_a = int(a_fields[1])
    rname_a = a_fields[2]
    pos_a = int(a_fields[3])
    mapq_a = int(a_fields[4])

    flag_b = int(b_fields[1])
    rname_b = b_fields[2]
    pos_b = int(b_fields[3])
    mapq_b = int(b_fields[4])

    is_unmapped_a = bool(flag_a & 0x4)
    is_unmapped_b = bool(flag_b & 0x4)
    is_duplicate_a = bool(flag_a & 0x400)
    is_duplicate_b = bool(flag_b & 0x400)

    if is_unmapped_a or is_unmapped_b:
        return False
    if is_duplicate_a or is_duplicate_b:
        return False
    if mapq_a == 0 or mapq_b == 0:
        return False
    if not passes_match_length(a_fields) or not passes_match_length(b_fields):
        return False
    if rname_a == rname_b:
        if check_distance:
            return abs(pos_a - pos_b) >= min_intra_dist
        return False
    return True

total_pairs = 0
informative_pairs = 0
chimeric_pairs_3d = 0
read_a = None

for line in sys.stdin:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6:
        continue
    qname = fields[0]

    if read_a is None:
        read_a = fields
        continue

    if qname == read_a[0]:
        total_pairs += 1
        if is_informative_pair(read_a, fields, check_distance=True):
            informative_pairs += 1
        if is_informative_pair(read_a, fields, check_distance=False):
            chimeric_pairs_3d += 1
        read_a = None
    else:
        read_a = fields

non_chimeric_pairs_3d = total_pairs - chimeric_pairs_3d
ratio_3d = float("nan") if total_pairs == 0 else chimeric_pairs_3d / total_pairs
informative_ratio = 0.0 if total_pairs == 0 else informative_pairs / total_pairs

print(f"Total pairs: {total_pairs}")
print(f"Informative pairs: {informative_pairs}")
print(f"Chimeric pairs (3D): {chimeric_pairs_3d}")
print(f"Non-chimeric pairs (3D): {non_chimeric_pairs_3d}")
print(f"3D ratio: {ratio_3d:.4f}")
print(f"Informative ratio: {informative_ratio:.4f}")

with open(f"{output_dir}/3d_ratio.txt", "w") as handle:
    handle.write(f"{ratio_3d:.4f}")
with open(f"{output_dir}/informative_pairs_ratio.txt", "w") as handle:
    handle.write(f"{informative_ratio:.4f}")
PY

    samtools view -f 1 -F 0x904 -q "$MAPQ" "$bam_path" | python3 "$metrics_script" "$output_dir" "$MIN_INTRA_DIST" "$MIN_MATCH_LEN"

    echo "[INFO] 3D ratio and informative pairs ratio calculation completed"
}

REFERENCE="output/assembly/final_assembly.fasta"
READS_1="output/readqc/hic/final_reads_1.fastq"
READS_2="output/readqc/hic/final_reads_2.fastq"
OUTPUT_DIR="output/alignment"
THREADS=80
BWA_OPTIONS="-5SP"
SAMTOOLS_FILTER="-F 0x900"
MAPQ=30
MIN_INTRA_DIST=10000
MIN_MATCH_LEN=30
SORT_MEMORY="1G"
TMP_DIR=""
KEEP_SAM=false
SKIP_METRICS=false
path=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p) path="$2"; shift 2 ;;
        -r|--reference) REFERENCE="$2"; shift 2 ;;
        -1|--reads1) READS_1="$2"; shift 2 ;;
        -2|--reads2) READS_2="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --bwa-options) BWA_OPTIONS="$2"; shift 2 ;;
        --samtools-filter) SAMTOOLS_FILTER="$2"; shift 2 ;;
        --mapq) MAPQ="$2"; shift 2 ;;
        --min-intra-dist) MIN_INTRA_DIST="$2"; shift 2 ;;
        --min-match-len) MIN_MATCH_LEN="$2"; shift 2 ;;
        --sort-memory) SORT_MEMORY="$2"; shift 2 ;;
        --tmp-dir) TMP_DIR="$2"; shift 2 ;;
        --keep-sam) KEEP_SAM=true; shift ;;
        --skip-metrics) SKIP_METRICS=true; shift ;;
        --print-defaults) print_defaults; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

for numeric in THREADS MAPQ MIN_INTRA_DIST MIN_MATCH_LEN; do
    value="${!numeric}"
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] ${numeric} must be a non-negative integer. Got: $value" >&2
        exit 1
    fi
done
if (( THREADS < 1 )); then
    echo "[ERROR] THREADS must be at least 1." >&2
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "[ERROR] Reference FASTA not found: $REFERENCE" >&2
    exit 1
fi
if [[ ! -f "$READS_1" ]]; then
    echo "[ERROR] Reads file not found: $READS_1" >&2
    exit 1
fi
if [[ ! -f "$READS_2" ]]; then
    echo "[ERROR] Reads file not found: $READS_2" >&2
    exit 1
fi

if [[ -z "$TMP_DIR" ]]; then
    if ! TMP_DIR="$(default_tmp_root)"; then
        echo "[ERROR] No writable temporary directory found. Set --tmp-dir, METAHICT_TMP_ROOT, or TMPDIR." >&2
        exit 1
    fi
fi
if [[ ! -d "$TMP_DIR" || ! -w "$TMP_DIR" ]]; then
    echo "[ERROR] Temporary directory is not writable: $TMP_DIR" >&2
    exit 1
fi

for tool in bwa samtools python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[ERROR] Required tool not found in PATH: $tool" >&2
        exit 1
    fi
done

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
echo "[INFO] Using reference: $REFERENCE"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/reference"

RUN_TMP="$(mktemp -d "${TMP_DIR%/}/metahict_alignment.XXXXXX")"
cleanup_tmp() {
    rm -rf "$RUN_TMP"
}
cleanup_then_return_status() {
    local status="$1"
    cleanup_tmp || true
    return "$status"
}
if [[ -n "${METAHICT_RESOURCE_FILE:-}" ]] && declare -F _metahict_resource_finish >/dev/null; then
    trap 'cleanup_then_return_status "$?"; _metahict_resource_finish' EXIT
else
    trap 'status="$?"; cleanup_tmp || true; exit "$status"' EXIT
fi

WORK_REFERENCE="${OUTPUT_DIR}/reference/$(basename "$REFERENCE")"
cp "$REFERENCE" "$WORK_REFERENCE"

read -r -a BWA_OPTIONS_ARRAY <<< "$BWA_OPTIONS"
read -r -a SAMTOOLS_FILTER_ARRAY <<< "$SAMTOOLS_FILTER"

MATCH_FILTER_SCRIPT="${RUN_TMP}/filter_by_match_length.py"
cat > "$MATCH_FILTER_SCRIPT" <<'PY'
import re
import sys

min_match_len = int(sys.argv[1])
cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")

def match_length(cigar):
    if cigar == "*":
        return 0
    return sum(int(n) for n, op in cigar_re.findall(cigar) if op in {"M", "=", "X"})

for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6:
        continue
    if min_match_len <= 0 or match_length(fields[5]) >= min_match_len:
        sys.stdout.write(line)
PY

echo "[INFO] Indexing copied reference with BWA: $WORK_REFERENCE"
bwa index "$WORK_REFERENCE"
if [[ ! -f "${WORK_REFERENCE}.bwt" ]]; then
    echo "[ERROR] BWA indexing failed: ${WORK_REFERENCE}.bwt not found" >&2
    exit 1
fi

echo "[INFO] Aligning reads with BWA MEM..."
echo "[INFO] Alignment filters: ${SAMTOOLS_FILTER}; MAPQ >= ${MAPQ}; aligned match length >= ${MIN_MATCH_LEN}"
if [[ "$KEEP_SAM" == true ]]; then
    bwa mem "${BWA_OPTIONS_ARRAY[@]}" -t "$THREADS" "$WORK_REFERENCE" "$READS_1" "$READS_2" > "$OUTPUT_DIR/map.sam"
    samtools view "${SAMTOOLS_FILTER_ARRAY[@]}" -h -q "$MAPQ" "$OUTPUT_DIR/map.sam" \
        | python3 "$MATCH_FILTER_SCRIPT" "$MIN_MATCH_LEN" \
        | samtools view -bS - \
        | samtools sort -n -@ "$THREADS" -m "$SORT_MEMORY" -T "${RUN_TMP}/sort" -o "$OUTPUT_DIR/sorted_map.bam" -
else
    bwa mem "${BWA_OPTIONS_ARRAY[@]}" -t "$THREADS" "$WORK_REFERENCE" "$READS_1" "$READS_2" \
        | samtools view "${SAMTOOLS_FILTER_ARRAY[@]}" -h -q "$MAPQ" - \
        | python3 "$MATCH_FILTER_SCRIPT" "$MIN_MATCH_LEN" \
        | samtools view -bS - \
        | samtools sort -n -@ "$THREADS" -m "$SORT_MEMORY" -T "${RUN_TMP}/sort" -o "$OUTPUT_DIR/sorted_map.bam" -
fi

if [[ ! -s "$OUTPUT_DIR/sorted_map.bam" ]]; then
    echo "[ERROR] Alignment failed: sorted BAM was not created." >&2
    exit 1
fi

echo "[INFO] Alignment completed successfully."

if [[ "$SKIP_METRICS" == false ]]; then
    echo "[INFO] Starting Hi-C quality metrics calculation..."
    calculate_hic_ratios "$OUTPUT_DIR/sorted_map.bam" "$OUTPUT_DIR"
else
    echo "[INFO] Skipping Hi-C quality metrics (--skip-metrics)."
fi

echo "[INFO] All alignment and analysis tasks completed successfully."
echo "[INFO] Results:"
echo "  - Sorted BAM file: $OUTPUT_DIR/sorted_map.bam"
if [[ "$KEEP_SAM" == true ]]; then
    echo "  - Intermediate SAM file: $OUTPUT_DIR/map.sam"
fi
if [[ "$SKIP_METRICS" == false ]]; then
    echo "  - 3D ratio: $OUTPUT_DIR/3d_ratio.txt"
    echo "  - Informative pairs ratio: $OUTPUT_DIR/informative_pairs_ratio.txt"
fi
