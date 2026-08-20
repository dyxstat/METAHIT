#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"
METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "$METAHICT_RESOURCE_OUTDIR" ]]; then
    metahict_resource_start "$METAHICT_RESOURCE_OUTDIR" "$(basename "$0")"
fi

set -euo pipefail
###########################################################################################################
# This script is a comprehensive solution to QC HiSeq reads in preparation for assembly and other
# operations. The main tasks are read trimming based on quality scores and producing FASTQC reports.
###########################################################################################################

# Function to display the help message
help_message () {
    echo ""
    echo "Usage: 1_preprocessing.sh [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
    echo "Options:"
    echo "    -p metahict path"
    echo "    -1 STR          forward fastq reads (.fastq or .fastq.gz)"
    echo "    -2 STR          reverse fastq reads (.fastq or .fastq.gz)"
    echo "    -o STR          output directory"
    echo "    --prefix STR    custom prefix for output files (default=base name extracted from input)"
    echo "    -t INT          number of threads (default=80)"
    echo "    --minlen INT    minimum length of reads after trimming (default=50)"
    echo "    --trimq INT     quality threshold for trimming (default=10)"
    echo "    --qtrim STR     quality trimming direction for BBDuk (default=r)"
    echo "    --ftl INT       trim bases from the left (default=10)"
    echo "    --xmx STR       memory for Java (default=80% of available memory)"
    echo "    --ftm INT       optional ftm value (default=5)"
    echo "    --ktrim STR     adapter trimming direction for BBDuk (default=r)"
    echo "    --k INT         adapter k-mer size for BBDuk (default=23)"
    echo "    --mink INT      minimum adapter k-mer size for BBDuk (default=11)"
    echo "    --hdist INT     Hamming distance for adapter k-mer matching (default=1)"
    echo "    --adapter-ref STR  adapter reference FASTA (default=BBTools packaged adapters.fa)"
    echo "    --dedup         perform deduplication with clumpify (default=false)"
    echo "    --skip-pre-qc-report   skip FastQC for input sequences"
    echo "    --skip-post-qc-report  skip FastQC for final sequences"
    echo "    --print-defaults       print default Module 1 parameter values and exit"
    echo ""
}

print_defaults () {
    cat <<EOF
Module 1 preprocessing defaults:
threads=80
prefix=derived from input R1 filename
dedup=false
minlen=50
trimq=10
qtrim=r
ftl=10
xmx=80% of available system memory
ftm=5
ktrim=r
k=23
mink=11
hdist=1
adapter_ref=BBTools packaged adapters.fa
skip_pre_qc_report=false
skip_post_qc_report=false
EOF
}

# Defaults
threads=80
minlen=50
trimq=10
qtrim=r
ftl=10
dedup=false
pre_qc_report=true
post_qc_report=true
out=""
reads_1=""
reads_2=""
prefix=""
path=""
ftm=5
ktrim=r
k=23
mink=11
hdist=1
adapter_ref=""

# Calculate default xmx as 80% of available memory
available_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
available_mem_gb=$((available_mem_kb / 1024 / 1024))

calculated_xmx=$((available_mem_gb * 80 / 100))
xmx="${calculated_xmx}g"

# Getopts
OPTS=$(getopt -o ht:o:1:2:p: --long help,print-defaults,skip-pre-qc-report,skip-post-qc-report,dedup,minlen:,trimq:,qtrim:,ftl:,xmx:,ftm:,ktrim:,k:,mink:,hdist:,adapter-ref:,prefix: -- "$@")
if [ $? -ne 0 ]; then help_message; exit 1; fi
eval set -- "$OPTS"

# Parse options
while true; do
    case "$1" in
        -p) path=$2; shift 2;;
        -t) threads=$2; shift 2;;
        -o) out=$2; shift 2;;
        -1) reads_1=$2; shift 2;;
        -2) reads_2=$2; shift 2;;
        --prefix) prefix=$2; shift 2;;
        --minlen) minlen=$2; shift 2;;
        --trimq) trimq=$2; shift 2;;
        --qtrim) qtrim=$2; shift 2;;
        --ftl) ftl=$2; shift 2;;
        --xmx) xmx=$2; shift 2;;
        --ftm) ftm=$2; shift 2;;
        --ktrim) ktrim=$2; shift 2;;
        --k) k=$2; shift 2;;
        --mink) mink=$2; shift 2;;
        --hdist) hdist=$2; shift 2;;
        --adapter-ref) adapter_ref=$2; shift 2;;
        --dedup) dedup=true; shift 1;;
        --skip-pre-qc-report) pre_qc_report=false; shift 1;;
        --skip-post-qc-report) post_qc_report=false; shift 1;;
        --print-defaults) print_defaults; exit 0;;
        -h|--help) help_message; exit 0;;
        --) shift; break;;
        *) echo "Invalid option: $1"; help_message; exit 1;;
    esac
done

# Use the pinned BBTools package supplied by the selected runtime environment.
# Prefer bbdukOld.sh when present because some BBTools releases route bbduk.sh
# to BBDukS, which can fail on large paired FASTQ streams.
bbduk_cmd="$(command -v bbdukOld.sh || command -v bbduk.sh || true)"
clumpify_cmd="$(command -v clumpify.sh || true)"
if [[ -z "$bbduk_cmd" || -z "$clumpify_cmd" ]]; then
    echo "Error: BBTools (bbduk.sh and clumpify.sh) is not available on PATH"
    exit 1
fi

# Set the adapter reference. A user-supplied path takes precedence; otherwise
# use the adapters file installed with the pinned Conda BBTools package.
if [ -n "$adapter_ref" ]; then
    ref="$adapter_ref"
else
    bbtools_prefix="$(cd "$(dirname "$bbduk_cmd")/.." && pwd)"
    ref="${bbtools_prefix}/share/bbmap/resources/adapters.fa"
    if [[ ! -f "$ref" ]]; then
        echo "Error: BBTools adapter reference not found: $ref"
        exit 1
    fi
fi

# Check if all required parameters are entered
if [ -z "$out" ] || [ -z "$reads_1" ] || [ -z "$reads_2" ] || [ -z "$path" ]; then
    help_message; exit 1
fi

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

# Check if output directory exists, if not, create it
mkdir -p "$out"

# Check if input read files exist in either .fastq or .fastq.gz formats
for read_file in "$reads_1" "$reads_2"; do
    if [[ ! -f "$read_file" ]] && [[ ! -f "${read_file}.gz" ]]; then
        echo "Error: $read_file or ${read_file}.gz not found"
        exit 1
    fi
done

# Generate base name if prefix is not set
if [ -z "$prefix" ]; then
    # Extract the base name without path and extension from first read file
    base_name1=$(basename "$reads_1")
    
    # Remove .fastq or .fastq.gz extension if present
    base_name1="${base_name1%.fastq}"
    base_name1="${base_name1%.fastq.gz}"
    
    # Remove the trailing _1 or _2 to get the base name
    base="${base_name1%_*}"
else
    # Use provided prefix
    base="$prefix"
fi

# FastQC report for input reads if requested
if [ "$pre_qc_report" = true ]; then
    mkdir -p "${out}/pre-QC_${base}_report"
    fastqc -q -t "$threads" -o "${out}/pre-QC_${base}_report" "$reads_1" "$reads_2"
    rm "${out}/pre-QC_${base}_report"/*zip
fi

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS

# Run BBDuk Steps with .gz output
echo "Running Adapter trimming with BBDuk"
"$bbduk_cmd" -Xmx$xmx in1="$reads_1" in2="$reads_2" out1="${out}/step1_adptrim_${base}_1.fastq.gz" out2="${out}/step1_adptrim_${base}_2.fastq.gz" \
    ref="$ref" ktrim="$ktrim" k="$k" mink="$mink" hdist="$hdist" minlen="$minlen" threads="$threads" tpe tbo

echo "Running Quality trimming with BBDuk"
"$bbduk_cmd" -Xmx$xmx in1="${out}/step1_adptrim_${base}_1.fastq.gz" in2="${out}/step1_adptrim_${base}_2.fastq.gz" \
    out1="${out}/step2_qualtrim_${base}_1.fastq.gz" out2="${out}/step2_qualtrim_${base}_2.fastq.gz" qtrim="$qtrim" trimq="$trimq" ftm="$ftm" \
    minlen="$minlen" threads="$threads"

echo "Trimming left bases with BBDuk"
"$bbduk_cmd" -Xmx$xmx in1="${out}/step2_qualtrim_${base}_1.fastq.gz" in2="${out}/step2_qualtrim_${base}_2.fastq.gz" \
    out1="${out}/step3_lefttrim_${base}_1.fastq.gz" out2="${out}/step3_lefttrim_${base}_2.fastq.gz" ftl="$ftl" threads="$threads"

# Deduplication step
if [ "$dedup" = true ]; then
    echo "Running deduplication with Clumpify"
    "$clumpify_cmd" in1="${out}/step3_lefttrim_${base}_1.fastq.gz" in2="${out}/step3_lefttrim_${base}_2.fastq.gz" \
        out1="${out}/final_${base}_1.fastq.gz" out2="${out}/final_${base}_2.fastq.gz" dedupe threads="$threads" 
else
    mv "${out}/step3_lefttrim_${base}_1.fastq.gz" "${out}/final_${base}_1.fastq.gz"
    mv "${out}/step3_lefttrim_${base}_2.fastq.gz" "${out}/final_${base}_2.fastq.gz"
fi

# FastQC report for final reads if requested
if [ "$post_qc_report" = true ]; then
    mkdir -p "${out}/post-QC_${base}_report"
    fastqc -t "$threads" -o "${out}/post-QC_${base}_report" "${out}/final_${base}_1.fastq.gz" "${out}/final_${base}_2.fastq.gz"
    rm "${out}/post-QC_${base}_report"/*zip
fi

echo "${base} preprocessing complete!"
