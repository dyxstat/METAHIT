#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"

# Parse arguments
if [ "$#" -lt 7 ]; then
    echo "Usage: $0 <normalization_method> -p <metahit_path> --bam <BAM file> --fasta <FASTA file> --out <output directory> --enzyme <enzyme name> [optional parameters]"
    exit 1
fi

COMMAND=$1
shift

BAM=""
FASTA=""
OUTDIR=""
ENZYME=""
path=""
METACC_MIN_SIGNAL=1
METACC_MIN_LEN=1000
METACC_MIN_MAPQ=30
METACC_MIN_MATCH=30

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p) path=$2; shift 2;;
        --bam) BAM="$2"; shift 2;;
        --fasta) FASTA="$2"; shift 2;;
        --out) OUTDIR="$2"; shift 2;;
        --enzyme) ENZYME="$2"; shift 2;;
        --metacc-min-signal) METACC_MIN_SIGNAL="$2"; shift 2;;
        --metacc-min-len) METACC_MIN_LEN="$2"; shift 2;;
        --metacc-min-mapq) METACC_MIN_MAPQ="$2"; shift 2;;
        --metacc-min-match) METACC_MIN_MATCH="$2"; shift 2;;
        *) echo "[ERROR] Unknown parameter: $1"; exit 1;;
    esac
done

if [ -z "$BAM" ] || [ -z "$FASTA" ] || [ -z "$OUTDIR" ] || [ -z "$ENZYME" ] || [ -z "$path" ]; then
    echo "[ERROR] Missing required arguments."
    exit 1
fi

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/tmp"

# Step 5a: Raw contact generation
python ${path}/5_contact/scripts/raw_contact.py "$BAM" "$ENZYME" "$FASTA" "$OUTDIR" "$METACC_MIN_MAPQ" "$METACC_MIN_LEN" "$METACC_MIN_MATCH" "$METACC_MIN_SIGNAL"
if [ $? -ne 0 ]; then
    echo "[ERROR] Raw contact matrix generation failed."
    exit 1
fi
echo "[INFO] Raw contact matrix generation completed successfully."

# Step 5b: Normalization
contig_file="${OUTDIR}/contig_info.csv"
contact_matrix_file="${OUTDIR}/Raw_contact_matrix.npz"
output_path="${OUTDIR}"

if [ "$COMMAND" == "metator" ] || [ "$COMMAND" == "hiczin" ]; then
    echo "[INFO] Merging coverage into contig info"
    coverage_file="${path}/m4/coverage_estimation/coverage.txt"
    merged_contig_file="${output_path}/contig_info_with_coverage.csv"

    python "${path}/5_contact/scripts/add_coverage.py" \
        --contig_info "$contig_file" \
        --coverage "$coverage_file" \
        --output "$merged_contig_file"

    contig_file="$merged_contig_file"
fi

python "${path}/5_contact/scripts/normalization.py" "$COMMAND" \
    --contig_file "$contig_file" \
    --contact_matrix_file "$contact_matrix_file" \
    --output_path "$output_path"

if [ $? -ne 0 ]; then
    echo "[ERROR] Normalization step '$COMMAND' failed."
    exit 1
fi

echo "[INFO] Normalization step '$COMMAND' completed successfully."
