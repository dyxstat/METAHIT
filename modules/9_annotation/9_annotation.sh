#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"

usage() {
    cat <<EOF
Usage: $0 -p PROJECT_PATH --bin BIN_DIR --outdir OUTDIR [options]

Required:
  -p, --project-path PATH          Path to METAHICT project root
  --bin PATH                       Directory containing genome bins
  --outdir PATH                    Output directory for annotation results

Options:
  -t, --threads INT                Number of CPU threads (default: 80)
  --pplacer-cpus INT               Number of CPUs for pplacer (default: same as --threads)
  --gtdbtk_db PATH                 Path to GTDB-Tk database (default: <project_path>/databases/gtdbtk_db/release220)
  --extension TEXT                 Genome file extension processed by GTDB-Tk (default: fa)
  --prefix TEXT                    Prefix for GTDB-Tk output files (default: gtdbtk)
  --skip-ani-screen                Skip ANI screening (default: true)
  --no-skip-ani-screen             Use GTDB-Tk ANI/Mash screening; requires --mash-db
  --mash-db PATH                   Mash database path used when ANI screening is enabled (default: empty)
  --min-perc-aa FLOAT              Minimum percentage of amino acids in the MSA (default: 10)
  --min-af FLOAT                   Minimum alignment fraction for species assignment (default: 0.5)
  --full-tree                      Use GTDB-Tk full bacterial tree (default: false)
  --scratch-dir PATH               Scratch directory for pplacer disk-backed memory reduction (default: empty)
  --tmp-dir PATH                   Temporary directory root (default: METAHICT_TMP_ROOT, TMPDIR, or /tmp)
  --force                          Continue if GTDB-Tk errors on a single genome (default: false)
  --keep-intermediates             Keep GTDB-Tk intermediate files (default: false)
  --debug                          Keep GTDB-Tk debug intermediates (default: false)
  --write-single-copy-genes        Output unaligned single-copy marker genes (default: false)
  --gtdbtk-extra-args TEXT         Additional native options passed to gtdbtk classify_wf, for example "--genes" or "--no_mash --mash_k 21" (default: empty)
  --print-defaults                 Print default Module 9 parameter values and exit
  -h, --help                       Show this help message and exit
EOF
}

print_defaults() {
    cat <<EOF
Module 9 annotation defaults:
threads=80
pplacer_cpus=same as threads
gtdbtk_db=<project_path>/databases/gtdbtk_db/release220
extension=fa
prefix=gtdbtk
skip_ani_screen=true
mash_db=
min_perc_aa=10
min_af=0.5
full_tree=false
scratch_dir=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
force=false
keep_intermediates=false
debug=false
write_single_copy_genes=false
gtdbtk_extra_args=
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

METAHICT_RESOURCE_OUTDIR="$(metahict_detect_outdir "$@" || true)"
if [[ -n "${METAHICT_RESOURCE_OUTDIR}" ]]; then
    metahict_resource_start "${METAHICT_RESOURCE_OUTDIR}" "$(basename "$0")"
fi

project_path=""
bin_dir=""
out_dir=""
threads=80
pplacer_cpus=""
gtdbtk_db=""
extension="fa"
prefix="gtdbtk"
skip_ani_screen=true
mash_db=""
min_perc_aa=10
min_af=0.5
full_tree=false
scratch_dir=""
tmp_dir="${METAHICT_TMP_ROOT:-${TMPDIR:-/tmp}}"
force=false
keep_intermediates=false
debug=false
write_single_copy_genes=false
gtdbtk_extra_args=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project-path)
            project_path="$2"
            shift 2
            ;;
        --bin)
            bin_dir="$2"
            shift 2
            ;;
        --outdir)
            out_dir="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        --pplacer-cpus)
            pplacer_cpus="$2"
            shift 2
            ;;
        --gtdbtk_db|--gtdbtk-db)
            gtdbtk_db="$2"
            shift 2
            ;;
        --extension)
            extension="$2"
            shift 2
            ;;
        --prefix)
            prefix="$2"
            shift 2
            ;;
        --skip-ani-screen)
            skip_ani_screen=true
            shift
            ;;
        --no-skip-ani-screen)
            skip_ani_screen=false
            shift
            ;;
        --mash-db)
            mash_db="$2"
            skip_ani_screen=false
            shift 2
            ;;
        --min-perc-aa)
            min_perc_aa="$2"
            shift 2
            ;;
        --min-af)
            min_af="$2"
            shift 2
            ;;
        --full-tree)
            full_tree=true
            shift
            ;;
        --scratch-dir)
            scratch_dir="$2"
            shift 2
            ;;
        --tmp-dir)
            tmp_dir="$2"
            shift 2
            ;;
        --force)
            force=true
            shift
            ;;
        --keep-intermediates)
            keep_intermediates=true
            shift
            ;;
        --debug)
            debug=true
            shift
            ;;
        --write-single-copy-genes)
            write_single_copy_genes=true
            shift
            ;;
        --gtdbtk-extra-args)
            gtdbtk_extra_args="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$project_path" || -z "$bin_dir" || -z "$out_dir" ]]; then
    echo "Error: Missing required argument(s)." >&2
    usage >&2
    exit 1
fi

if [[ -z "$pplacer_cpus" ]]; then
    pplacer_cpus="$threads"
fi

if [[ -z "$gtdbtk_db" ]]; then
    gtdbtk_db="${project_path}/databases/gtdbtk_db/release220"
fi

if [[ "$skip_ani_screen" == false && -z "$mash_db" ]]; then
    echo "Error: --mash-db is required when --no-skip-ani-screen is used." >&2
    exit 1
fi

mkdir -p "$out_dir" "$tmp_dir"
export GTDBTK_DATA_PATH="$gtdbtk_db"
export TMPDIR="$tmp_dir"

if command -v gtdbtk >/dev/null 2>&1; then
    gtdbtk_bin="$(command -v gtdbtk)"
elif [[ -x "${project_path}/conda_envs/gtdbtk-2.4.0/bin/gtdbtk" ]]; then
    gtdbtk_bin="${project_path}/conda_envs/gtdbtk-2.4.0/bin/gtdbtk"
else
    echo "Error: gtdbtk executable not found. Install GTDB-Tk or provide a METAHICT project with conda_envs/gtdbtk-2.4.0." >&2
    exit 1
fi
export PATH="$(dirname "$gtdbtk_bin"):${PATH}"

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
echo "[INFO] Using GTDB-Tk executable: $gtdbtk_bin"
echo "[INFO] Using GTDB-Tk database: $GTDBTK_DATA_PATH"
echo "[INFO] Running GTDB-Tk classify_wf with $threads threads and $pplacer_cpus pplacer CPUs..."

cmd=(
    "$gtdbtk_bin" classify_wf
    --genome_dir "$bin_dir"
    --out_dir "$out_dir"
    --extension "$extension"
    --prefix "$prefix"
    --cpus "$threads"
    --pplacer_cpus "$pplacer_cpus"
    --min_perc_aa "$min_perc_aa"
    --min_af "$min_af"
    --tmpdir "$tmp_dir"
)

if [[ "$skip_ani_screen" == true ]]; then
    cmd+=(--skip_ani_screen)
else
    cmd+=(--mash_db "$mash_db")
fi

if [[ "$full_tree" == true ]]; then
    cmd+=(--full_tree)
fi
if [[ -n "$scratch_dir" ]]; then
    mkdir -p "$scratch_dir"
    cmd+=(--scratch_dir "$scratch_dir")
fi
if [[ "$force" == true ]]; then
    cmd+=(--force)
fi
if [[ "$keep_intermediates" == true ]]; then
    cmd+=(--keep_intermediates)
fi
if [[ "$debug" == true ]]; then
    cmd+=(--debug)
fi
if [[ "$write_single_copy_genes" == true ]]; then
    cmd+=(--write_single_copy_genes)
fi
if [[ -n "$gtdbtk_extra_args" ]]; then
    read -r -a extra_args <<< "$gtdbtk_extra_args"
    cmd+=("${extra_args[@]}")
fi

printf '[INFO] Command:'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

echo "Annotation (GTDB-Tk) step completed successfully."
