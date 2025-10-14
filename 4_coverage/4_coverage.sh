#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
# Coverage Estimation Script using BBMap and jgi_summarize_bam_contig_depths.

set -e
set -o pipefail

# Function to print informational messages
function echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

# Function to print error messages
function echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

# Display usage if not enough arguments
if [ "$#" -lt 8 ]; then
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir [-m Java_heap_memory]"
    exit 1
fi

# Set default Java heap memory to 80% of available memory
default_memory=$(( $(free -m | awk '/^Mem:/{print $2}') * 80 / 100 ))
JAVA_HEAP="-Xmx${default_memory}m"

threads=20

# Read in the arguments
while getopts "p:1:2:r:o:m:t:" opt; do
    case $opt in
        p) path=$OPTARG ;;
        1) reads_1=$OPTARG ;;
        2) reads_2=$OPTARG ;;
        r) ref=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        m) JAVA_HEAP="-Xmx$OPTARG" ;;
        t) threads=$OPTARG ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

# Check if required parameters are set
if [[ -z "$reads_1" || -z "$reads_2" || -z "$ref" || -z "$output_dir" ]]; then
    echo_error "Missing required parameters"
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir [-m Java_heap_memory]"
    exit 1
fi

# Clear and recreate output directory
if [[ -d "$output_dir" ]]; then
    echo_info "Removing existing output directory: $output_dir"
    rm -rf "$output_dir"
fi
mkdir -p "$output_dir"

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS

# Run BBMap with memory setting and capture output for debugging
echo_info "Running BBMap..."
${path}/external/bbmap/bbmap.sh \
    in1="$reads_1" in2="$reads_2" ref="$ref" \
    out="$output_dir/SG_map.sam" bamscript="$output_dir/bs.sh" \
    threads="$threads" $JAVA_HEAP > "$output_dir/bbmap.log" 2>&1

# Check if BBMap generated the SAM file and bamscript successfully
if [[ ! -f "$output_dir/SG_map.sam" || ! -f "$output_dir/bs.sh" ]]; then
    echo_error "BBMap failed to generate SAM or bamscript. Check $output_dir/bbmap.log for details."
    exit 1
fi

# Run SAM to sorted BAM conversion using the generated script
echo_info "Converting SAM to sorted BAM with bs.sh..."
chmod +x "$output_dir/bs.sh"
sh "$output_dir/bs.sh" > "$output_dir/bs.log" 2>&1

# Check if sorted BAM was successfully generated
if [[ ! -f "$output_dir/SG_map_sorted.bam" ]]; then
    echo_error "Sorted BAM (SG_map_sorted.bam) not found. Check $output_dir/bs.log for errors."
    exit 1
fi

# Run jgi_summarize_bam_contig_depths and log output
echo_info "Running jgi_summarize_bam_contig_depths..."
${path}/4_coverage/scripts/jgi_summarize_bam_contig_depths \
    --outputDepth "$output_dir/coverage.txt" \
    --pairedContigs "$output_dir/pair.txt" \
    "$output_dir/SG_map_sorted.bam" > "$output_dir/jgi_summarize.log" 2>&1

# Check if coverage.txt was created
if [[ -f "$output_dir/coverage.txt" ]]; then
    echo_info "Coverage.txt generated successfully at $output_dir/coverage.txt"
else
    echo_error "Failed to generate coverage.txt. Check $output_dir/jgi_summarize.log for details."
    exit 1
fi

