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
Usage: $0 -p PROJECT_PATH --fasta BIN_FASTA --enzyme ENZYMES --outdir OUTDIR --hic1 HIC_R1 --hic2 HIC_R2 [options]

Required:
  -p, --project-path PATH          Path to METAHICT project root or modules directory
  --fasta PATH                     Input bin FASTA
  --enzyme TEXT                    Restriction enzyme names, comma-separated when multiple enzymes are used
  --outdir PATH                    Output directory
  --hic1 PATH                      Forward Hi-C reads
  --hic2 PATH                      Reverse Hi-C reads

Options:
  --bam PATH                       Optional Hi-C BAM. If omitted, reads are aligned internally.
  -t, --threads INT                Number of CPU threads (default: 80)
  -m, --memory TEXT                Memory limit passed to the module (default: 80% of available memory)
  -r, --resolution INT             Segment length for scaffold visualization (default: 10000)
  --min-contig-len INT             Minimum contig length retained for scaffolding (default: 5000)
  --bwa-options TEXT               BWA MEM options for internal alignments (default: -5SP)
  --samtools-filter TEXT           samtools view filter for internal alignments (default: -F 0x900)
  --sort-memory TEXT               Memory per samtools sort thread for internal alignments (default: 1G)
  --metacc-min-mapq INT            Minimum MAPQ for MetaCC contact generation (default: 30)
  --metacc-min-len INT             Minimum contig length for MetaCC contact generation (default: 1000)
  --metacc-min-match INT           Minimum aligned match length for MetaCC contact generation (default: 30)
  --metacc-min-signal INT          Minimum contact signal for MetaCC contact generation (default: 2)
  --bin3c-min-mapq INT             Minimum MAPQ for bin3C-compatible contact generation (default: 60)
  --bin3c-min-len INT              Minimum contig length for bin3C-compatible contact generation (default: 1000)
  --bin3c-min-match INT            Minimum aligned match length for bin3C-compatible contact generation (default: 10)
  --bin3c-min-signal INT           Minimum contact signal for bin3C-compatible contact generation (default: 5)
  --yahs-resolutions TEXT          YaHS scaffolding resolution list; empty uses YaHS automatic selection (default: empty)
  --yahs-min-mapq INT              YaHS minimum mapping quality (default: 10)
  --yahs-min-contig-len INT        YaHS minimum contig length to scaffold (default: 0)
  --yahs-rounds INT                YaHS rounds at each resolution level (default: 1)
  --yahs-no-contig-ec              Disable YaHS contig error correction (default: false)
  --yahs-no-scaffold-ec            Disable YaHS scaffold error correction (default: false)
  --yahs-no-mem-check              Disable YaHS runtime memory check (default: false)
  --yahs-extra-args TEXT           Additional options passed directly to YaHS (default: empty)
  --normcc-thres FLOAT             NormCC denoising threshold (default: 0.05)
  --heatmap-max-image INT          Maximum heatmap image dimension before downsampling (default: 5000)
  --skip-checkm2                   Skip CheckM2 quality evaluation (default: false)
  --checkm2_db PATH                Path to custom CheckM2 database
  --tmp-dir PATH                   Temporary directory root (default: METAHICT_TMP_ROOT, TMPDIR, or /tmp)
  --keep-temp                      Keep temporary files for debugging (default: false)
  --print-defaults                 Print default Module 8 parameter values and exit
  -h, --help                       Show this help message and exit
EOF
}

print_defaults() {
    cat <<EOF
Module 8 scaffolding defaults:
threads=80
memory=80% of available system memory
resolution=10000
min_contig_len=5000
bwa_options=-5SP
samtools_filter=-F 0x900
sort_memory=1G
metacc_min_mapq=30
metacc_min_len=1000
metacc_min_match=30
metacc_min_signal=2
bin3c_min_mapq=60
bin3c_min_len=1000
bin3c_min_match=10
bin3c_min_signal=5
yahs_resolutions=
yahs_min_mapq=10
yahs_min_contig_len=0
yahs_rounds=1
yahs_no_contig_ec=false
yahs_no_scaffold_ec=false
yahs_no_mem_check=false
yahs_extra_args=
normcc_thres=0.05
heatmap_max_image=5000
skip_checkm2=false
checkm2_db=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
keep_temp=false
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

available_mem_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)
available_mem_gb=$((available_mem_kb / 1024 / 1024))
calculated_mem=$((available_mem_gb * 80 / 100))

project_path=""
fasta=""
bam=""
enzyme=""
outdir=""
hic1=""
hic2=""
threads=80
memory="${calculated_mem}g"
resolution=10000
min_contig_len=5000
bwa_options="-5SP"
samtools_filter="-F 0x900"
sort_memory="1G"
metacc_min_mapq=30
metacc_min_len=1000
metacc_min_match=30
metacc_min_signal=2
bin3c_min_mapq=60
bin3c_min_len=1000
bin3c_min_match=10
bin3c_min_signal=5
yahs_resolutions=""
yahs_min_mapq=10
yahs_min_contig_len=0
yahs_rounds=1
yahs_no_contig_ec=0
yahs_no_scaffold_ec=0
yahs_no_mem_check=0
yahs_extra_args=""
normcc_thres=0.05
heatmap_max_image=5000
skip_checkm2=0
checkm2_db=""
tmp_dir="${METAHICT_TMP_ROOT:-${TMPDIR:-/tmp}}"
keep_temp=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project-path)
            project_path="$2"; shift 2;;
        --fasta)
            fasta="$2"; shift 2;;
        --bam)
            bam="$2"; shift 2;;
        --enzyme)
            enzyme="$2"; shift 2;;
        --outdir)
            outdir="$2"; shift 2;;
        --hic1)
            hic1="$2"; shift 2;;
        --hic2)
            hic2="$2"; shift 2;;
        -t|--threads)
            threads="$2"; shift 2;;
        -m|--memory)
            memory="$2"; shift 2;;
        -r|--resolution)
            resolution="$2"; shift 2;;
        --min-contig-len)
            min_contig_len="$2"; shift 2;;
        --bwa-options)
            bwa_options="$2"; shift 2;;
        --samtools-filter)
            samtools_filter="$2"; shift 2;;
        --sort-memory)
            sort_memory="$2"; shift 2;;
        --metacc-min-mapq)
            metacc_min_mapq="$2"; shift 2;;
        --metacc-min-len)
            metacc_min_len="$2"; shift 2;;
        --metacc-min-match)
            metacc_min_match="$2"; shift 2;;
        --metacc-min-signal)
            metacc_min_signal="$2"; shift 2;;
        --bin3c-min-mapq)
            bin3c_min_mapq="$2"; shift 2;;
        --bin3c-min-len)
            bin3c_min_len="$2"; shift 2;;
        --bin3c-min-match)
            bin3c_min_match="$2"; shift 2;;
        --bin3c-min-signal)
            bin3c_min_signal="$2"; shift 2;;
        --yahs-resolutions)
            yahs_resolutions="$2"; shift 2;;
        --yahs-min-mapq)
            yahs_min_mapq="$2"; shift 2;;
        --yahs-min-contig-len)
            yahs_min_contig_len="$2"; shift 2;;
        --yahs-rounds)
            yahs_rounds="$2"; shift 2;;
        --yahs-no-contig-ec)
            yahs_no_contig_ec=1; shift;;
        --yahs-no-scaffold-ec)
            yahs_no_scaffold_ec=1; shift;;
        --yahs-no-mem-check)
            yahs_no_mem_check=1; shift;;
        --yahs-extra-args)
            yahs_extra_args="$2"; shift 2;;
        --normcc-thres)
            normcc_thres="$2"; shift 2;;
        --heatmap-max-image)
            heatmap_max_image="$2"; shift 2;;
        --skip-checkm2)
            skip_checkm2=1; shift;;
        --checkm2_db)
            checkm2_db="$2"; shift 2;;
        --tmp-dir)
            tmp_dir="$2"; shift 2;;
        --keep-temp)
            keep_temp=1; shift;;
        --print-defaults)
            print_defaults
            exit 0;;
        -h|--help)
            usage
            exit 0;;
        *)
            echo_error "Unknown parameter passed: $1"
            usage
            exit 1;;
    esac
done

if [[ -z "${project_path}" || -z "${fasta}" || -z "${enzyme}" || -z "${outdir}" || -z "${hic1}" || -z "${hic2}" ]]; then
    echo_error "Missing one or more required parameters."
    usage
    exit 1
fi

for input in "${fasta}" "${hic1}" "${hic2}"; do
    if [[ ! -e "${input}" ]]; then
        echo_error "Missing input: ${input}"
        exit 1
    fi
done
if [[ -n "${bam}" && ! -e "${bam}" ]]; then
    echo_error "Missing BAM input: ${bam}"
    exit 1
fi

if [[ -d "${project_path}/8_scaffolding/scripts" ]]; then
    module_root="${project_path}"
elif [[ -d "${project_path}/modules/8_scaffolding/scripts" ]]; then
    module_root="${project_path}/modules"
else
    echo_error "Could not locate modules/8_scaffolding/scripts from -p ${project_path}"
    exit 1
fi

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: ${free_mem}"

mkdir -p "${outdir}" "${tmp_dir}"
export MPLCONFIGDIR="${outdir}/.matplotlib"
mkdir -p "${MPLCONFIGDIR}"
scaffolding_script="${module_root}/8_scaffolding/scripts/scaffolding.py"

cmd=(
    python "${scaffolding_script}"
    -p "${module_root}"
    --fasta "${fasta}"
    --enzyme "${enzyme}"
    --outdir "${outdir}"
    --hic1 "${hic1}"
    --hic2 "${hic2}"
    -t "${threads}"
    -m "${memory}"
    -r "${resolution}"
    --min-contig-len "${min_contig_len}"
    "--bwa-options=${bwa_options}"
    "--samtools-filter=${samtools_filter}"
    --sort-memory "${sort_memory}"
    --metacc-min-mapq "${metacc_min_mapq}"
    --metacc-min-len "${metacc_min_len}"
    --metacc-min-match "${metacc_min_match}"
    --metacc-min-signal "${metacc_min_signal}"
    --bin3c-min-mapq "${bin3c_min_mapq}"
    --bin3c-min-len "${bin3c_min_len}"
    --bin3c-min-match "${bin3c_min_match}"
    --bin3c-min-signal "${bin3c_min_signal}"
    --yahs-min-mapq "${yahs_min_mapq}"
    --yahs-min-contig-len "${yahs_min_contig_len}"
    --yahs-rounds "${yahs_rounds}"
    --normcc-thres "${normcc_thres}"
    --heatmap-max-image "${heatmap_max_image}"
    --tmp-dir "${tmp_dir}"
)

if [[ -n "${bam}" ]]; then
    cmd+=(--bam "${bam}")
fi
if [[ -n "${yahs_resolutions}" ]]; then
    cmd+=(--yahs-resolutions "${yahs_resolutions}")
fi
if [[ "${yahs_no_contig_ec}" -eq 1 ]]; then
    cmd+=(--yahs-no-contig-ec)
fi
if [[ "${yahs_no_scaffold_ec}" -eq 1 ]]; then
    cmd+=(--yahs-no-scaffold-ec)
fi
if [[ "${yahs_no_mem_check}" -eq 1 ]]; then
    cmd+=(--yahs-no-mem-check)
fi
if [[ -n "${yahs_extra_args}" ]]; then
    cmd+=(--yahs-extra-args "${yahs_extra_args}")
fi
if [[ "${skip_checkm2}" -eq 1 ]]; then
    cmd+=(--skip-checkm2)
fi
if [[ -n "${checkm2_db}" ]]; then
    cmd+=(--checkm2_db "${checkm2_db}")
fi
if [[ "${keep_temp}" -eq 1 ]]; then
    cmd+=(--keep-temp)
fi

echo_info "Running scaffolding.py"
"${cmd[@]}"
echo_info "Scaffolding step completed successfully."
