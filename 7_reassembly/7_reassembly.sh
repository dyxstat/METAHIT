#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

# Script to perform reassembly for Metahit
if [ "$#" -lt 18 ]; then
    echo "Usage: $0 -p <METAHIT_PATH> --bin <BIN> --assembly <ASSEMBLY> --hic1 <HIC1> --hic2 <HIC2> --sg1 <SG1> --sg2 <SG2> --bam <BAM> --outdir <OUTDIR> -t <THREADS> -m <MEMORY>"
    exit 1
fi

# Parse arguments properly
while [[ $# -gt 0 ]]; do
    case $1 in
        -p) PATH_DIR="$2"; shift 2 ;;
        --bin) BIN="$2"; shift 2 ;;
        --assembly) ASSEMBLY="$2"; shift 2 ;;
        --hic1) HIC1="$2"; shift 2 ;;
        --hic2) HIC2="$2"; shift 2 ;;
        --sg1) SG1="$2"; shift 2 ;;
        --sg2) SG2="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -m) MEMORY="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "$PATH_DIR" || -z "$BIN" || -z "$ASSEMBLY" || -z "$HIC1" || -z "$HIC2" || -z "$SG1" || -z "$SG2" || -z "$BAM" || -z "$OUTDIR" || -z "$THREADS" || -z "$MEMORY" ]]; then
    echo "Error: Missing required arguments"
    exit 1
fi

# Create contig sort list
mkdir -p "$OUTDIR"
CONTIG_SORT_LIST="${OUTDIR}/contig_sort.txt"
echo "Creating contig sort list: $CONTIG_SORT_LIST"
if [ ! -f "${ASSEMBLY}.fai" ]; then
    samtools faidx "$ASSEMBLY"
fi
sort -k2,2nr "${ASSEMBLY}.fai" | cut -f1 > "$CONTIG_SORT_LIST"

# Path to the reassembly Python script
REASSEMBLY_SCRIPT="${PATH_DIR}/7_reassembly/scripts/7_reassembly.py"

# Run the reassembly script with provided arguments
python "$REASSEMBLY_SCRIPT" \
    --bin "$BIN" \
    --hic1 "$HIC1" \
    --hic2 "$HIC2" \
    --sg1 "$SG1" \
    --sg2 "$SG2" \
    --bam "$BAM" \
    --outdir "$OUTDIR" \
    -p "$PATH_DIR" \
    -c "$CONTIG_SORT_LIST" \
    -t "$THREADS" \
    -m "$MEMORY"

# Check if the reassembly process succeeded
if [ $? -ne 0 ]; then
    echo "Error: Reassembly process failed."
    exit 1
fi

echo "Reassembly completed successfully."