#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"
METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "$METAHICT_RESOURCE_OUTDIR" ]]; then
    metahict_resource_start "$METAHICT_RESOURCE_OUTDIR" "$(basename "$0")"
fi

set -euo pipefail

default_memory_gb() {
    awk '/^MemTotal:/ {printf "%d\n", ($2 / 1024 / 1024) * 0.8}' /proc/meminfo
}

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
    cat <<EOF
Module 2 assembly defaults:
threads=80
memory=80% of available system memory
min_len=1000
assembler=megahit
megahit_k_min=21
megahit_k_max=141
megahit_k_step=12
megahit_merge_level=20,0.95
metaspades_k_list=21,33,55
flye_method=--nano-raw
tmp_dir=${tmp_dir:-METAHICT_TMP_ROOT, TMPDIR, or /tmp}
skip_quast=false
keep_temp=false
EOF
}

help_message() {
    cat <<'EOF'

Usage: metahict assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir
Options:
    -p STR              METAHICT project path or modules path
    -1 STR              Forward FASTQ reads
    -2 STR              Reverse FASTQ reads
    -o STR              Output directory
    -t INT              Number of CPU threads (default=80)
    -m INT              Memory in GB (default=80% of available system memory)
    -l INT              Minimum assembled contig length (default=1000)
    -h, --help          Show this help message
    --print-defaults    Print default parameter values and exit

Assembler selection:
    --megahit           Assemble with MEGAHIT (default)
    --metaspades        Assemble with metaSPAdes
    --metaflye          Assemble with metaFlye

MEGAHIT options:
    --k-min INT         Minimum k-mer size (default=21)
    --k-max INT         Maximum k-mer size (default=141)
    --k-step INT        K-mer step size (default=12)
    --merge-level STR   Merge level (default=20,0.95)

metaSPAdes options:
    --k-list STR        Comma-separated k-mer sizes (default=21,33,55)

metaFlye options:
    --flye-method STR   Flye read mode (default=--nano-raw)

Runtime/output options:
    --tmp-dir STR       Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)
    --skip-quast        Skip QUAST quality-control report
    --keep-temp         Keep temporary assembler files

EOF
}

mem="$(default_memory_gb)"
threads=80
out=""
reads_1=""
reads_2=""
path=""
min_len=1000
merge_level="20,0.95"
flye_method="--nano-raw"
k_min=21
k_max=141
k_step=12
k_list="21,33,55"
tmp_dir=""
skip_quast=false
keep_temp=false
assembler="megahit"
assembler_flags=0

OPTS=$(getopt -o hp:t:m:o:1:2:l: --long help,print-defaults,megahit,metaspades,metaflye,flye-method:,k-min:,k-max:,k-step:,k-list:,merge-level:,tmp-dir:,skip-quast,keep-temp -- "$@")
if [[ $? != 0 ]]; then help_message; exit 1; fi
eval set -- "$OPTS"

while true; do
    case "$1" in
        -h|--help) help_message; exit 0 ;;
        --print-defaults) print_defaults; exit 0 ;;
        -p) path="$2"; shift 2 ;;
        -t) threads="$2"; shift 2 ;;
        -m) mem="$2"; shift 2 ;;
        -o) out="$2"; shift 2 ;;
        -1) reads_1="$2"; shift 2 ;;
        -2) reads_2="$2"; shift 2 ;;
        -l) min_len="$2"; shift 2 ;;
        --megahit) assembler="megahit"; assembler_flags=$((assembler_flags + 1)); shift ;;
        --metaspades) assembler="metaspades"; assembler_flags=$((assembler_flags + 1)); shift ;;
        --metaflye) assembler="metaflye"; assembler_flags=$((assembler_flags + 1)); shift ;;
        --flye-method) flye_method="$2"; shift 2 ;;
        --k-min) k_min="$2"; shift 2 ;;
        --k-max) k_max="$2"; shift 2 ;;
        --k-step) k_step="$2"; shift 2 ;;
        --k-list) k_list="$2"; shift 2 ;;
        --merge-level) merge_level="$2"; shift 2 ;;
        --tmp-dir) tmp_dir="$2"; shift 2 ;;
        --skip-quast) skip_quast=true; shift ;;
        --keep-temp) keep_temp=true; shift ;;
        --) shift; break ;;
        *) echo "Unknown option: $1"; help_message; exit 1 ;;
    esac
done

if (( assembler_flags > 1 )); then
    echo "[ERROR] Choose only one assembler: --megahit, --metaspades, or --metaflye." >&2
    exit 1
fi

if [[ -z "$path" || -z "$out" || -z "$reads_1" || -z "$reads_2" ]]; then
    help_message
    exit 1
fi

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

if ! [[ "$threads" =~ ^[0-9]+$ ]] || (( threads < 1 )); then
    echo "[ERROR] -t/--threads must be a positive integer." >&2
    exit 1
fi
if ! [[ "$mem" =~ ^[0-9]+$ ]] || (( mem < 1 )); then
    echo "[ERROR] -m/--memory must be a positive integer in GB." >&2
    exit 1
fi
if ! [[ "$min_len" =~ ^[0-9]+$ ]] || (( min_len < 1 )); then
    echo "[ERROR] -l/--min-len must be a positive integer." >&2
    exit 1
fi

if [[ -z "$tmp_dir" ]]; then
    if ! tmp_dir="$(default_tmp_root)"; then
        echo "[ERROR] No writable temporary directory found. Set --tmp-dir, METAHICT_TMP_ROOT, or TMPDIR." >&2
        exit 1
    fi
fi
if [[ ! -d "$tmp_dir" || ! -w "$tmp_dir" ]]; then
    echo "[ERROR] Temporary directory is not writable: $tmp_dir" >&2
    exit 1
fi

if [[ -d "${path}/2_assembly/scripts" ]]; then
    MODULE_ROOT="$path"
elif [[ -d "${path}/modules/2_assembly/scripts" ]]; then
    MODULE_ROOT="${path}/modules"
else
    echo "[ERROR] Cannot locate module scripts from -p path: $path" >&2
    exit 1
fi
SOFT="${MODULE_ROOT}/2_assembly/scripts"

echo "[INFO] Checking dependencies..."
if [[ "$assembler" == "metaspades" ]] && ! command -v metaspades.py &> /dev/null; then
    echo "[ERROR] metaSPAdes is not installed or not in the PATH." >&2
    exit 1
fi
if [[ "$assembler" == "megahit" ]] && ! command -v megahit &> /dev/null; then
    echo "[ERROR] MEGAHIT is not installed or not in the PATH." >&2
    exit 1
fi
if [[ "$assembler" == "metaflye" ]] && ! command -v flye &> /dev/null; then
    echo "[ERROR] Flye is not installed or not in the PATH." >&2
    exit 1
fi
if [[ "$skip_quast" == false ]] && ! command -v quast.py &> /dev/null; then
    echo "[ERROR] QUAST is not installed or not in the PATH." >&2
    exit 1
fi
if ! command -v samtools &> /dev/null; then
    echo "[ERROR] samtools is not installed or not in the PATH." >&2
    exit 1
fi
if [[ ! -f "${SOFT}/rm_short_contigs.py" ]]; then
    echo "[ERROR] Script 'rm_short_contigs.py' not found in ${SOFT}." >&2
    exit 1
fi

mkdir -p "$out" || { echo "[ERROR] Cannot create output directory: $out" >&2; exit 1; }
mkdir -p "${out}/.matplotlib"
export TMPDIR="$tmp_dir"
export MPLCONFIGDIR="${out}/.matplotlib"

if [[ "$assembler" == "metaspades" ]]; then
    IFS=',' read -ra K_ARRAY <<< "$k_list"
    for k in "${K_ARRAY[@]}"; do
        if ! [[ $k =~ ^[0-9]+$ ]]; then
            echo "[ERROR] metaSPAdes k-mer size '$k' is not an integer." >&2
            exit 1
        fi
        if (( k % 2 == 0 )); then
            echo "[ERROR] metaSPAdes k-mer size '$k' is not odd." >&2
            exit 1
        fi
        if (( k >= 128 )); then
            echo "[ERROR] metaSPAdes k-mer size '$k' must be less than 128." >&2
            exit 1
        fi
    done
fi

if [[ "$assembler" == "megahit" ]]; then
    if ! [[ $k_min =~ ^[0-9]+$ ]] || ! [[ $k_max =~ ^[0-9]+$ ]] || ! [[ $k_step =~ ^[0-9]+$ ]]; then
        echo "[ERROR] MEGAHIT --k-min, --k-max, and --k-step must be integers." >&2
        exit 1
    fi
    if (( k_min % 2 == 0 || k_max % 2 == 0 )); then
        echo "[ERROR] MEGAHIT --k-min and --k-max must be odd." >&2
        exit 1
    fi
    if (( k_step % 2 != 0 )); then
        echo "[ERROR] MEGAHIT --k-step must be even." >&2
        exit 1
    fi
    if (( k_min < 15 || k_max > 255 || k_min > k_max )); then
        echo "[ERROR] Invalid MEGAHIT k-mer range: k-min=$k_min, k-max=$k_max." >&2
        exit 1
    fi
fi

run_tmp=""
cleanup_tmp() {
    if [[ -n "$run_tmp" && "$keep_temp" == false ]]; then
        rm -rf "$run_tmp"
    fi
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

if [[ "$assembler" == "metaspades" ]]; then
    echo "ASSEMBLING WITH metaSPAdes"
    run_tmp="$(mktemp -d "${tmp_dir%/}/metahict_metaspades.XXXXXX")"
    metaspades.py --tmp-dir "$run_tmp" -k "$k_list" -t "$threads" -m "$mem" -o "${out}/metaspades" -1 "$reads_1" -2 "$reads_2"
    if [[ ! -f "${out}/metaspades/scaffolds.fasta" ]]; then
        echo "[ERROR] metaSPAdes assembly failed." >&2
        exit 1
    fi
    "${SOFT}/rm_short_contigs.py" "$min_len" "${out}/metaspades/scaffolds.fasta" > "${out}/final_assembly.fasta"
fi

if [[ "$assembler" == "megahit" ]]; then
    echo "ASSEMBLING WITH MEGAHIT"
    run_tmp="$(mktemp -d "${tmp_dir%/}/metahict_megahit.XXXXXX")"
    mem_bytes=$((mem * 1024 * 1024 * 1024))
    megahit -1 "$reads_1" -2 "$reads_2" -o "${out}/megahit" --min-contig-len "$min_len" --k-min "$k_min" --k-max "$k_max" --k-step "$k_step" --merge-level "$merge_level" -t "$threads" -m "$mem_bytes" --tmp-dir "$run_tmp"
    if [[ ! -f "${out}/megahit/final.contigs.fa" ]]; then
        echo "[ERROR] MEGAHIT assembly failed." >&2
        exit 1
    fi
    cp "${out}/megahit/final.contigs.fa" "${out}/final_assembly.fasta"
fi

if [[ "$assembler" == "metaflye" ]]; then
    echo "[INFO] Assembling with metaFlye..."
    run_tmp="$(mktemp -d "${tmp_dir%/}/metahict_metaflye.XXXXXX")"
    flye "$flye_method" "$reads_1" "$reads_2" --meta -t "$threads" --out-dir "${out}/metaflye"
    if [[ ! -f "${out}/metaflye/assembly.fasta" ]]; then
        echo "[ERROR] Flye assembly failed." >&2
        exit 1
    fi
    cp "${out}/metaflye/assembly.fasta" "${out}/final_assembly.fasta"
fi

if [[ ! -s "${out}/final_assembly.fasta" ]]; then
    echo "[ERROR] Final assembly failed." >&2
    exit 1
fi

if [[ "$skip_quast" == false ]]; then
    echo "RUNNING ASSEMBLY QC WITH QUAST"
    quast.py -t "$threads" -o "${out}/QUAST_out" "${out}/final_assembly.fasta"
else
    echo "[INFO] Skipping QUAST QC (--skip-quast)."
fi
samtools faidx "${out}/final_assembly.fasta"
