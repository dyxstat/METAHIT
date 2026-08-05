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

usage() {
    cat <<EOF
Usage: $0 <normalization_method> -p PROJECT_PATH --bam BAM --fasta FASTA --out OUTDIR --enzyme ENZYMES [options]

Normalization methods:
  raw, normcc, hiczin, bin3c, metator

Required:
  -p, --project-path PATH          Path to METAHICT project root or modules directory
  --bam PATH                       Queryname-sorted Hi-C BAM file
  --fasta PATH                     Assembly FASTA
  --out, --outdir PATH             Output directory
  --enzyme TEXT                    Restriction enzyme names, comma-separated when multiple enzymes are used

Options:
  --metacc-min-signal INT          Minimum off-diagonal contact signal for retained contigs (default: 1)
  --metacc-min-len INT             Minimum contig length for raw contact generation (default: 1000)
  --metacc-min-mapq INT            Minimum MAPQ for raw contact generation (default: 30)
  --metacc-min-match INT           Minimum terminal aligned match length for raw contact generation (default: 30)
  --spurious-contact-percent FLOAT Percentile cutoff for spurious-contact denoising (default: 5)
  --thres FLOAT                    Alias for --spurious-contact-percent
  --coverage-file PATH             Coverage file from Module 4; required for hiczin and metator (default: empty)
  --epsilon FLOAT                  Epsilon used by hiczin, bin3c, and metator normalization (default: 1)
  --max-iter INT                   Maximum iterations for bin3c matrix balancing (default: 1000)
  --tol FLOAT                      Convergence tolerance for bin3c matrix balancing (default: 1e-6)
  --print-defaults                 Print default Module 5 parameter values and exit
  -h, --help                       Show this help message and exit
EOF
}

print_defaults() {
    cat <<EOF
Module 5 contact defaults:
normalization_method=required
metacc_min_signal=1
metacc_min_len=1000
metacc_min_mapq=30
metacc_min_match=30
spurious_contact_percent=5
coverage_file=
epsilon=1
max_iter=1000
tol=1e-6
EOF
}

if [[ $# -gt 0 ]]; then
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --print-defaults)
            print_defaults
            exit 0
            ;;
    esac
fi

if [[ $# -lt 1 ]]; then
    echo_error "Missing normalization method"
    usage
    exit 1
fi

METHOD="$1"
shift

case "${METHOD}" in
    raw|normcc|hiczin|bin3c|metator)
        ;;
    *)
        echo_error "Unsupported normalization method: ${METHOD}"
        usage
        exit 1
        ;;
esac

BAM=""
FASTA=""
OUTDIR=""
ENZYME=""
PROJECT_PATH=""
METACC_MIN_SIGNAL=1
METACC_MIN_LEN=1000
METACC_MIN_MAPQ=30
METACC_MIN_MATCH=30
SPURIOUS_CONTACT_PERCENT=5
COVERAGE_FILE=""
EPSILON=1
MAX_ITER=1000
TOL=1e-6

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project-path)
            PROJECT_PATH="$2"
            shift 2
            ;;
        --bam)
            BAM="$2"
            shift 2
            ;;
        --fasta)
            FASTA="$2"
            shift 2
            ;;
        --out|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --enzyme)
            ENZYME="$2"
            shift 2
            ;;
        --metacc-min-signal)
            METACC_MIN_SIGNAL="$2"
            shift 2
            ;;
        --metacc-min-len)
            METACC_MIN_LEN="$2"
            shift 2
            ;;
        --metacc-min-mapq)
            METACC_MIN_MAPQ="$2"
            shift 2
            ;;
        --metacc-min-match)
            METACC_MIN_MATCH="$2"
            shift 2
            ;;
        --spurious-contact-percent|--thres)
            SPURIOUS_CONTACT_PERCENT="$2"
            shift 2
            ;;
        --coverage-file)
            COVERAGE_FILE="$2"
            shift 2
            ;;
        --epsilon)
            EPSILON="$2"
            shift 2
            ;;
        --max-iter)
            MAX_ITER="$2"
            shift 2
            ;;
        --tol)
            TOL="$2"
            shift 2
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
            echo_error "Unknown parameter: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${BAM}" || -z "${FASTA}" || -z "${OUTDIR}" || -z "${ENZYME}" || -z "${PROJECT_PATH}" ]]; then
    echo_error "Missing required arguments."
    usage
    exit 1
fi

for input in "${BAM}" "${FASTA}"; do
    if [[ ! -e "${input}" ]]; then
        echo_error "Missing input: ${input}"
        exit 1
    fi
done

if [[ "${METHOD}" == "hiczin" || "${METHOD}" == "metator" ]]; then
    if [[ -z "${COVERAGE_FILE}" ]]; then
        echo_error "--coverage-file is required when method is ${METHOD}"
        exit 1
    fi
    if [[ ! -e "${COVERAGE_FILE}" ]]; then
        echo_error "Missing coverage file: ${COVERAGE_FILE}"
        exit 1
    fi
fi

if [[ -d "${PROJECT_PATH}/5_contact/scripts" ]]; then
    MODULE_ROOT="${PROJECT_PATH}"
elif [[ -d "${PROJECT_PATH}/modules/5_contact/scripts" ]]; then
    MODULE_ROOT="${PROJECT_PATH}/modules"
else
    echo_error "Could not locate modules/5_contact/scripts from -p ${PROJECT_PATH}"
    exit 1
fi

METAHICT_ROOT="$(dirname "${MODULE_ROOT}")"
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: ${free_mem}"

if [[ -d "${OUTDIR}" ]]; then
    echo_info "Removing existing output directory contents: ${OUTDIR}"
    find "${OUTDIR}" -mindepth 1 \
        ! -name "module.log" \
        ! -name "resources.txt" \
        ! -name ".resources_peak.*" \
        -exec rm -rf {} +
fi
mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/tmp"

echo_info "Generating raw contact matrix..."
python "${MODULE_ROOT}/5_contact/scripts/raw_contact.py" \
    "${BAM}" \
    "${ENZYME}" \
    "${FASTA}" \
    "${OUTDIR}" \
    "${METACC_MIN_MAPQ}" \
    "${METACC_MIN_LEN}" \
    "${METACC_MIN_MATCH}" \
    "${METACC_MIN_SIGNAL}"
echo_info "Raw contact matrix generation completed successfully."

CONTIG_FILE="${OUTDIR}/contig_info.csv"
CONTACT_MATRIX_FILE="${OUTDIR}/Raw_contact_matrix.npz"
OUTPUT_PATH="${OUTDIR}"

if [[ "${METHOD}" == "hiczin" || "${METHOD}" == "metator" ]]; then
    echo_info "Merging coverage into contig info"
    MERGED_CONTIG_FILE="${OUTPUT_PATH}/contig_info_with_coverage.csv"
    python "${MODULE_ROOT}/5_contact/scripts/add_coverage.py" \
        --contig_info "${CONTIG_FILE}" \
        --coverage "${COVERAGE_FILE}" \
        --output "${MERGED_CONTIG_FILE}"
    CONTIG_FILE="${MERGED_CONTIG_FILE}"
fi

NORMALIZATION_CMD=(
    python "${MODULE_ROOT}/5_contact/scripts/normalization.py" "${METHOD}"
    --contig_file "${CONTIG_FILE}"
    --contact_matrix_file "${CONTACT_MATRIX_FILE}"
    --output_path "${OUTPUT_PATH}"
    --thres "${SPURIOUS_CONTACT_PERCENT}"
    --min_len "${METACC_MIN_LEN}"
    --min_signal "${METACC_MIN_SIGNAL}"
)

case "${METHOD}" in
    hiczin|metator)
        NORMALIZATION_CMD+=(--epsilon "${EPSILON}")
        ;;
    bin3c)
        NORMALIZATION_CMD+=(--epsilon "${EPSILON}" --max_iter "${MAX_ITER}" --tol "${TOL}")
        ;;
esac

echo_info "Running ${METHOD} normalization..."
"${NORMALIZATION_CMD[@]}"
echo_info "Normalization step '${METHOD}' completed successfully."
