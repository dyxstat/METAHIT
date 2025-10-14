#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
set -e
set -o pipefail

function echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}
function echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

# Defaults
available_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
available_mem_gb=$((available_mem_kb / 1024 / 1024))
calculated_mem=$((available_mem_gb * 80 / 100))
mem="${calculated_mem}g"
threads=40
resolution=1000

# If insufficient arguments, show help
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 -p metahit_path --fasta bin.fa [--bam hic_mapped.bam] --enzyme DpnII --outdir scaffolding_output --hic1 hic_R1.fq --hic2 hic_R2.fq -t 40 -m 200g -r 10000"
    exit 1
fi

# Read arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p) path=$2; shift 2;;
        --fasta) fasta=$2; shift 2;;
        --bam) bam=$2; shift 2;;
        --enzyme) enzyme=$2; shift 2;;
        --outdir) outdir=$2; shift 2;;
        --hic1) hic1=$2; shift 2;;
        --hic2) hic2=$2; shift 2;;
        -t) threads=$2; shift 2;;
        -m) mem=$2; shift 2;;
        -r) resolution=$2; shift 2;;
        *) echo_error "Unknown parameter passed: $1"; exit 1;;
    esac
done

# Check required inputs (BAM is now optional)
if [[ -z "$path" || -z "$fasta" || -z "$enzyme" || -z "$outdir" || -z "$hic1" || -z "$hic2" ]]; then
    echo_error "Missing one or more required parameters."
    echo "Usage: $0 -p metahit_path --fasta bin.fa [--bam hic_mapped.bam] --enzyme DpnII --outdir output_dir --hic1 hic_R1.fq --hic2 hic_R2.fq [-t threads] [-m mem] [-r resolution]"
    exit 1
fi

# Check script existence
scaffolding_script="${path}/8_scaffolding/scripts/scaffolding.py"
if [ ! -f "$scaffolding_script" ]; then
    echo_error "scaffolding.py not found at $scaffolding_script"
    exit 1
fi

# Build command with optional BAM
echo_info "Running scaffolding.py"
cmd="python \"$scaffolding_script\" -p \"$path\" --fasta \"$fasta\" --enzyme \"$enzyme\" --outdir \"$outdir\" --hic1 \"$hic1\" --hic2 \"$hic2\" -t \"$threads\" -m \"$mem\" -r \"$resolution\""

# Add BAM only if provided
if [[ -n "$bam" ]]; then
    cmd="$cmd --bam \"$bam\""
fi

eval $cmd

if [ $? -eq 0 ]; then
    echo_info "Scaffolding step completed successfully."
else
    echo_error "Scaffolding step failed."
    exit 1
fi