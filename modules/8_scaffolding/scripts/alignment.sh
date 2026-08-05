#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 -r REFERENCE -1 READS_1 -2 READS_2 -o OUTPUT_DIR [options]

Options:
  -p PATH                         METAHICT modules path
  -r, --reference PATH            Reference FASTA
  -1, --reads1 PATH               Forward reads
  -2, --reads2 PATH               Reverse reads
  -o, --output PATH               Output directory
  -t, --threads INT               Number of CPU threads (default: 80)
  --bwa-options TEXT              BWA MEM options (default: -5SP)
  --samtools-filter TEXT          samtools view filter options (default: -F 0x900)
  --sort-memory TEXT              Memory per samtools sort thread (default: 1G)
  -h, --help                      Show this help message and exit
EOF
}

REFERENCE=""
READS_1=""
READS_2=""
OUTPUT_DIR=""
THREADS=80
BWA_OPTIONS="-5SP"
SAMTOOLS_FILTER="-F 0x900"
SORT_MEMORY="1G"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            shift 2;;
        -r|--reference)
            REFERENCE="$2"; shift 2;;
        -1|--reads1)
            READS_1="$2"; shift 2;;
        -2|--reads2)
            READS_2="$2"; shift 2;;
        -o|--output)
            OUTPUT_DIR="$2"; shift 2;;
        -t|--threads)
            THREADS="$2"; shift 2;;
        --bwa-options)
            BWA_OPTIONS="$2"; shift 2;;
        --samtools-filter)
            SAMTOOLS_FILTER="$2"; shift 2;;
        --sort-memory)
            SORT_MEMORY="$2"; shift 2;;
        -h|--help)
            usage
            exit 0;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage
            exit 1;;
    esac
done

if [[ -z "${REFERENCE}" || -z "${READS_1}" || -z "${READS_2}" || -z "${OUTPUT_DIR}" ]]; then
    echo "[ERROR] Missing required alignment arguments." >&2
    usage
    exit 1
fi

for input in "${REFERENCE}" "${READS_1}" "${READS_2}"; do
    if [[ ! -e "${input}" ]]; then
        echo "[ERROR] Missing input: ${input}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

echo "[INFO] Using reference: ${REFERENCE}"
echo "[INFO] Indexing reference genome with BWA..."
bwa index "${REFERENCE}"

if [[ ! -f "${REFERENCE}.bwt" ]]; then
    echo "[ERROR] BWA indexing failed: ${REFERENCE}.bwt not found" >&2
    exit 1
fi

echo "[INFO] Aligning reads with BWA MEM..."
bwa mem ${BWA_OPTIONS} -t "${THREADS}" "${REFERENCE}" "${READS_1}" "${READS_2}" > "${OUTPUT_DIR}/map.sam"

echo "[INFO] Converting SAM to BAM..."
samtools view ${SAMTOOLS_FILTER} -bS "${OUTPUT_DIR}/map.sam" > "${OUTPUT_DIR}/unsorted_map.bam"

echo "[INFO] Sorting BAM by read name..."
samtools sort -n -@ "${THREADS}" -m "${SORT_MEMORY}" "${OUTPUT_DIR}/unsorted_map.bam" -o "${OUTPUT_DIR}/sorted_map.bam"
rm -f "${OUTPUT_DIR}/unsorted_map.bam"
echo "[INFO] Alignment completed successfully!"
