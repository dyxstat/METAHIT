#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
set -e
set -o pipefail
set -x

# Default values for options
SAMTOOLS_FILTER="-F 0x900"
THREADS=20
# add BWA and samtools memory parameters if available
# Help message
usage() {
    echo "   -p metahit path"
    echo ""
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -r, --reference     Path to the reference genome (default: output/assembly/final_assembly.fasta)"
    echo "  -1, --reads1        Path to the first reads file (default: output/readqc/hic/final_reads_1.fastq)"
    echo "  -2, --reads2        Path to the second reads file (default: output/readqc/hic/final_reads_2.fastq)"
    echo "  -o, --output        Output directory (default: output/alignment)"
    echo "  -t, --threads       Number of threads to use (default: 20)"
    echo "  --samtools-filter   samtools view filter options (default: '0x904')"
    echo "  -h, --help          Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 --reference ref.fa --reads1 reads_1.fq --reads2 reads_2.fq --output alignment_output --threads 20 --samtools-filter '-F 0x904'"
    echo ""
}

# Parse command-line arguments
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -p) path=$2; shift 2;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -1|--reads1)
            READS_1="$2"
            shift 2
            ;;
        -2|--reads2)
            READS_2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --samtools-filter)
            SAMTOOLS_FILTER="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --) # end of all options
            shift
            break
            ;;
        -*|--*)
            echo "Unknown option $1"
            usage
            exit 1
            ;;
        *) # positional argument
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Set default values if not set
REFERENCE=${REFERENCE:-"output/assembly/final_assembly.fasta"}
READS_1=${READS_1:-"output/readqc/hic/final_reads_1.fastq"}
READS_2=${READS_2:-"output/readqc/hic/final_reads_2.fastq"}
OUTPUT_DIR=${OUTPUT_DIR:-"output/alignment"}

if [ ! -f "$REFERENCE" ]; then
    echo "[ERROR] Reference FASTA not found: $REFERENCE"
    exit 1
fi

echo "[INFO] Using reference: $REFERENCE"

# Define tool paths
BWA_PATH="bwa"
SAMTOOLS_PATH="samtools"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Index the reference genome
echo "Indexing reference genome with BWA..."
$BWA_PATH index "$REFERENCE"

if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "[ERROR] BWA indexing failed: ${REFERENCE}.bwt not found"
    exit 1
fi

# Align reads with BWA MEM
echo "Aligning reads with BWA MEM..."
$BWA_PATH mem -5SP -t "$THREADS" "$REFERENCE" "$READS_1" "$READS_2" > "$OUTPUT_DIR/map.sam"

# Convert SAM to BAM
echo "Converting SAM to BAM..."
$SAMTOOLS_PATH view $SAMTOOLS_FILTER -bS "$OUTPUT_DIR/map.sam" > "$OUTPUT_DIR/unsorted_map.bam"

# Sort BAM by read name
echo "Sorting BAM by read name..."
$SAMTOOLS_PATH sort -n "$OUTPUT_DIR/unsorted_map.bam" -o "$OUTPUT_DIR/sorted_map.bam"
rm "$OUTPUT_DIR/unsorted_map.bam"
echo "Alignment completed successfully!"
