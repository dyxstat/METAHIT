#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
set -e
set -o pipefail
set -x

# Default values for options
SAMTOOLS_FILTER="-F 0x900"
THREADS=20

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
    echo "  --samtools-filter   samtools view filter options (default: '0x900')"
    echo "  -h, --help          Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 --reference ref.fa --reads1 reads_1.fq --reads2 reads_2.fq --output alignment_output --threads 20 --samtools-filter '-F 0x900'"
    echo ""
}

# Combined function to calculate both 3D ratio and informative pairs ratio
calculate_hic_ratios() {
    local bam_path="$1"
    local output_dir="$2"
    
    echo "[INFO] Calculating 3D ratio and informative pairs ratio..."
    echo "[INFO] Extracting and processing primary paired alignments..."
    
    # Use samtools filtering and Python to calculate both ratios in one pass
    samtools view -f 1 -F 0x904 -q 30 "$bam_path" | python3 -c "
import sys

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
    
    is_informative = True
    if is_unmapped_a or is_unmapped_b:
        is_informative = False
    elif is_duplicate_a or is_duplicate_b:
        is_informative = False
    elif mapq_a == 0 or mapq_b == 0:
        is_informative = False
    elif rname_a == rname_b:
        if check_distance:
            if abs(pos_a - pos_b) < 10000:
                is_informative = False
        else:
            is_informative = False
    
    return is_informative

total_pairs = 0
informative_pairs = 0
chimeric_pairs_3d = 0
read_a = None

for line in sys.stdin:
    fields = line.strip().split('\t')
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

print(f'Total pairs: {total_pairs}')
print(f'Informative pairs: {informative_pairs}')
print(f'Chimeric pairs (3D): {chimeric_pairs_3d}')
print(f'Non-chimeric pairs (3D): {total_pairs - chimeric_pairs_3d}')

if total_pairs - chimeric_pairs_3d > 0:
    ratio_3d = chimeric_pairs_3d / (total_pairs - chimeric_pairs_3d)
else:
    ratio_3d = float('inf')

if total_pairs > 0:
    informative_ratio = informative_pairs / total_pairs
else:
    informative_ratio = 0.0

print(f'3D ratio: {ratio_3d:.4f}')
print(f'Informative ratio: {informative_ratio:.4f}')

# Save to files
with open(\"$output_dir/3d_ratio.txt\", 'w') as f:
    if ratio_3d == float('inf'):
        f.write('inf')
    else:
        f.write(f'{ratio_3d:.4f}')

with open(\"$output_dir/informative_pairs_ratio.txt\", 'w') as f:
    f.write(f'{informative_ratio:.4f}')
"
    
    echo "[INFO] 3D ratio and informative pairs ratio calculation completed"
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

# Calculate Hi-C quality metrics using the sorted BAM file
echo "[INFO] Starting Hi-C quality metrics calculation..."
calculate_hic_ratios "$OUTPUT_DIR/sorted_map.bam" "$OUTPUT_DIR"

echo "[INFO] All alignment and analysis tasks completed successfully!"
echo "[INFO] Results:"
echo "  - Sorted BAM file: $OUTPUT_DIR/sorted_map.bam"
echo "  - 3D ratio: $OUTPUT_DIR/3d_ratio.txt"
echo "  - Informative pairs ratio: $OUTPUT_DIR/informative_pairs_ratio.txt"