#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"
METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "${METAHICT_RESOURCE_OUTDIR}" ]]; then
    metahict_resource_start "${METAHICT_RESOURCE_OUTDIR}" "$(basename "$0")"
fi

echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

available_memory_mb() {
    free -m | awk '/^Mem:/{print $2}'
}

usage() {
    cat <<EOF
Usage: $0 -p PROJECT_PATH -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir [options]

Required:
  -p, --project-path PATH       Path to METAHICT project root or modules directory
  -1 PATH                       Forward shotgun reads
  -2 PATH                       Reverse shotgun reads
  -r, --reference PATH          Assembly FASTA
  -o, --output PATH             Output directory

Options:
  -t, --threads INT             Number of CPU threads (default: 80)
  -m, --memory VALUE            Java heap for BBMap, for example 64g or 64000m (default: 80% of available system memory)
  --percent-identity FLOAT      Minimum end-to-end percent identity for coverage counting (default: 97)
  --min-mapq INT                Minimum mapping quality for coverage counting (default: 0)
  --weight-mapq FLOAT           Weight per-base depth by mapping quality (default: 0.0; disabled)
  --include-edge-bases          Include read-length edge bases in depth and variance calculation (default: false)
  --max-edge-bases INT          Maximum edge length excluded when edge bases are not included (default: 75)
  --min-contig-length INT       Minimum contig length included by JGI coverage summarization (default: 0; disabled)
  --min-contig-depth FLOAT      Minimum contig depth used by JGI coverage summarization (default: 0; disabled)
  --bbmap-extra-args STRING     Additional options passed directly to bbmap.sh (default: empty)
  --tmp-dir PATH                Temporary directory root (default: METAHICT_TMP_ROOT, TMPDIR, or /tmp)
  --keep-sam                    Keep intermediate SG_map.sam (default: false)
  --keep-temp                   Keep temporary working/index files (default: false)
  --print-defaults              Print default Module 4 parameter values and exit
  -h, --help                    Show this help message and exit
EOF
}

print_defaults() {
    cat <<EOF
Module 4 coverage defaults:
threads=80
memory=80% of available system memory
percent_identity=97
min_mapq=0
weight_mapq=0.0
include_edge_bases=false
max_edge_bases=75
min_contig_length=0
min_contig_depth=0
bbmap_extra_args=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
keep_sam=false
keep_temp=false
EOF
}

THREADS=80
MEMORY=""
PERCENT_IDENTITY=97
MIN_MAPQ=0
WEIGHT_MAPQ=0.0
INCLUDE_EDGE_BASES=0
MAX_EDGE_BASES=75
MIN_CONTIG_LENGTH=0
MIN_CONTIG_DEPTH=0
BBMAP_EXTRA_ARGS=""
TMP_ROOT="${METAHICT_TMP_ROOT:-${TMPDIR:-/tmp}}"
KEEP_SAM=0
KEEP_TEMP=0

PROJECT_PATH=""
READS_1=""
READS_2=""
REF=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project-path)
            PROJECT_PATH="$2"
            shift 2
            ;;
        -1)
            READS_1="$2"
            shift 2
            ;;
        -2)
            READS_2="$2"
            shift 2
            ;;
        -r|--reference)
            REF="$2"
            shift 2
            ;;
        -o|--output|--outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        --percent-identity)
            PERCENT_IDENTITY="$2"
            shift 2
            ;;
        --min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;
        --weight-mapq)
            WEIGHT_MAPQ="$2"
            shift 2
            ;;
        --include-edge-bases)
            INCLUDE_EDGE_BASES=1
            shift
            ;;
        --max-edge-bases)
            MAX_EDGE_BASES="$2"
            shift 2
            ;;
        --min-contig-length)
            MIN_CONTIG_LENGTH="$2"
            shift 2
            ;;
        --min-contig-depth)
            MIN_CONTIG_DEPTH="$2"
            shift 2
            ;;
        --bbmap-extra-args)
            BBMAP_EXTRA_ARGS="$2"
            shift 2
            ;;
        --tmp-dir)
            TMP_ROOT="$2"
            shift 2
            ;;
        --keep-sam)
            KEEP_SAM=1
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=1
            shift
            ;;
        --print-defaults)
            print_defaults
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${PROJECT_PATH}" || -z "${READS_1}" || -z "${READS_2}" || -z "${REF}" || -z "${OUTPUT_DIR}" ]]; then
    echo_error "Missing required parameters"
    usage
    exit 1
fi

for input in "${READS_1}" "${READS_2}" "${REF}"; do
    if [[ ! -e "${input}" ]]; then
        echo_error "Missing input: ${input}"
        exit 1
    fi
done

if [[ -z "${MEMORY}" ]]; then
    MEMORY="$(( $(available_memory_mb) * 80 / 100 ))m"
fi
JAVA_HEAP="-Xmx${MEMORY}"

if [[ -x "${PROJECT_PATH}/4_coverage/scripts/jgi_summarize_bam_contig_depths" ]]; then
    MODULE_ROOT="${PROJECT_PATH}"
elif [[ -x "${PROJECT_PATH}/modules/4_coverage/scripts/jgi_summarize_bam_contig_depths" ]]; then
    MODULE_ROOT="${PROJECT_PATH}/modules"
else
    echo_error "Could not locate modules/4_coverage/scripts/jgi_summarize_bam_contig_depths from -p ${PROJECT_PATH}"
    exit 1
fi

JGI_SUMMARIZE="${MODULE_ROOT}/4_coverage/scripts/jgi_summarize_bam_contig_depths"
BBMAP="$(command -v bbmap.sh || true)"
if [[ -z "${BBMAP}" ]]; then
    echo_error "Could not locate bbmap.sh on PATH. Install the version-pinned BBMap dependency."
    exit 1
fi

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: ${free_mem}"

if [[ -d "${OUTPUT_DIR}" ]]; then
    echo_info "Removing existing output directory contents: ${OUTPUT_DIR}"
    find "${OUTPUT_DIR}" -mindepth 1 \
        ! -name "module.log" \
        ! -name "resources.txt" \
        ! -name ".resources_peak.*" \
        -exec rm -rf {} +
fi
mkdir -p "${OUTPUT_DIR}"

mkdir -p "${TMP_ROOT}"
RUN_TMP="$(mktemp -d "${TMP_ROOT%/}/metahict_coverage.XXXXXX")"
BBMAP_INDEX_DIR="${OUTPUT_DIR}/bbmap_index"
mkdir -p "${BBMAP_INDEX_DIR}"

cleanup() {
    if [[ "${KEEP_TEMP}" != "1" ]]; then
        rm -rf "${RUN_TMP}" "${BBMAP_INDEX_DIR}"
        rm -f "${OUTPUT_DIR}/bs.sh"
    fi
    if [[ "${KEEP_SAM}" != "1" ]]; then
        rm -f "${OUTPUT_DIR}/SG_map.sam"
    fi
}

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS
export TMPDIR="${RUN_TMP}"

BBMAP_CMD=(
    "${BBMAP}"
    "in1=${READS_1}"
    "in2=${READS_2}"
    "ref=${REF}"
    "path=${BBMAP_INDEX_DIR}"
    "out=${OUTPUT_DIR}/SG_map.sam"
    "bamscript=${OUTPUT_DIR}/bs.sh"
    "threads=${THREADS}"
    "${JAVA_HEAP}"
)
if [[ -n "${BBMAP_EXTRA_ARGS}" ]]; then
    read -r -a BBMAP_EXTRA_ARRAY <<< "${BBMAP_EXTRA_ARGS}"
    BBMAP_CMD+=("${BBMAP_EXTRA_ARRAY[@]}")
fi

echo_info "Running BBMap..."
"${BBMAP_CMD[@]}" 2>&1 | tee "${OUTPUT_DIR}/bbmap.log"

if [[ ! -f "${OUTPUT_DIR}/SG_map.sam" || ! -f "${OUTPUT_DIR}/bs.sh" ]]; then
    echo_error "BBMap failed to generate SAM or bamscript. Check ${OUTPUT_DIR}/bbmap.log for details."
    exit 1
fi

echo_info "Converting SAM to sorted BAM with bs.sh..."
chmod +x "${OUTPUT_DIR}/bs.sh"
sh "${OUTPUT_DIR}/bs.sh" 2>&1 | tee "${OUTPUT_DIR}/bs.log"

if [[ ! -f "${OUTPUT_DIR}/SG_map_sorted.bam" ]]; then
    echo_error "Sorted BAM (SG_map_sorted.bam) not found. Check ${OUTPUT_DIR}/bs.log for errors."
    exit 1
fi

JGI_CMD=(
    "${JGI_SUMMARIZE}"
    --outputDepth "${OUTPUT_DIR}/coverage.txt"
    --pairedContigs "${OUTPUT_DIR}/pair.txt"
    --percentIdentity "${PERCENT_IDENTITY}"
    --minMapQual "${MIN_MAPQ}"
    --weightMapQual "${WEIGHT_MAPQ}"
    --maxEdgeBases "${MAX_EDGE_BASES}"
)
if [[ "${INCLUDE_EDGE_BASES}" == "1" ]]; then
    JGI_CMD+=(--includeEdgeBases)
fi
if awk -v value="${MIN_CONTIG_LENGTH}" 'BEGIN {exit !(value > 0)}'; then
    JGI_CMD+=(--minContigLength "${MIN_CONTIG_LENGTH}")
fi
if awk -v value="${MIN_CONTIG_DEPTH}" 'BEGIN {exit !(value > 0)}'; then
    JGI_CMD+=(--minContigDepth "${MIN_CONTIG_DEPTH}")
fi
JGI_CMD+=("${OUTPUT_DIR}/SG_map_sorted.bam")

echo_info "Running jgi_summarize_bam_contig_depths..."
"${JGI_CMD[@]}" 2>&1 | tee "${OUTPUT_DIR}/jgi_summarize.log"

if [[ -f "${OUTPUT_DIR}/coverage.txt" ]]; then
    echo_info "Coverage.txt generated successfully at ${OUTPUT_DIR}/coverage.txt"
else
    echo_error "Failed to generate coverage.txt. Check ${OUTPUT_DIR}/jgi_summarize.log for details."
    exit 1
fi

cleanup
