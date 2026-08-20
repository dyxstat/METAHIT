#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../1_preprocessing/resource_logging.sh"

usage() {
    cat <<'EOF'
Usage: 10_MGE.sh -p PROJECT_PATH --combined COMBINED_FASTA --contact CONTACT_NPZ --raw-contact RAW_CONTACT_NPZ --outdir OUTDIR [options]

Required:
  -p, --metahict-path PATH          Path to METAHICT installation directory
  --combined PATH                  Combined contigs FASTA with bin/unmapped prefixes
  --contact PATH                   Normalized contact matrix (.npz)
  --raw-contact PATH               Raw contact matrix (.npz) used for raw Hi-C support
  --outdir PATH                    Output directory

Options:
  -t, --threads INT                Number of CPU threads (default: 80)
  --genomad-db PATH                geNomad database path (default: <project_path>/databases/genomad_db/genomad_db)
  --genomad-splits INT            geNomad MMseqs2 split count (default: 8)
  --genomad-sensitivity FLOAT      geNomad MMseqs2 sensitivity (default: 4.2)
  --genomad-cleanup                Delete geNomad intermediate files (default: true)
  --no-genomad-cleanup             Keep geNomad intermediate files
  --genomad-restart                Overwrite existing geNomad intermediate files (default: false)
  --genomad-preset TEXT            geNomad filtering preset: default, conservative, or relaxed (default: default)
  --genomad-min-score FLOAT        geNomad minimum virus/plasmid score when preset is default (default: 0.7)
  --genomad-max-fdr FLOAT          geNomad maximum FDR when preset is default (default: 0.1)
  --genomad-extra-args TEXT        Additional native options passed to genomad end-to-end (default: empty)
  --association-filter TEXT        Association filter: zscore, fixed, percentage, or raw-support-only (default: zscore)
  --zscore-threshold FLOAT         Minimum Z-score for filtered MGE-MAG associations (default: 0.5)
  --fixed-contact-threshold FLOAT  Minimum contact strength for fixed association filtering (default: 0)
  --top-percent FLOAT              Top percent of association contact strengths retained by percentage filtering (default: 50)
  --min-raw-contacts FLOAT         Minimum raw Hi-C read-pair support for candidate associations (default: 2)
  --ccfind-terminal-fragment-size INT  ccfind terminal fragment size (default: 500)
  --ccfind-min-identity INT            ccfind minimum terminal-alignment percent identity (default: 94)
  --ccfind-min-aligned-length INT      ccfind minimum terminal-alignment length (default: 50)
  --min-contact-strength FLOAT     Minimum positive contact strength counted for MGE-MAG associations (default: 0)
  --tmp-dir PATH                   Temporary directory root (default: METAHICT_TMP_ROOT, TMPDIR, or /tmp)
  --print-defaults                 Print default Module 10 parameter values and exit
  -h, --help                       Show this help message and exit
EOF
}

print_defaults() {
    cat <<'EOF'
Module 10 MGE defaults:
threads=80
genomad_db=<project_path>/databases/genomad_db/genomad_db
genomad_splits=8
genomad_sensitivity=4.2
genomad_cleanup=true
genomad_restart=false
genomad_preset=default
genomad_min_score=0.7
genomad_max_fdr=0.1
genomad_extra_args=
association_filter=zscore
zscore_threshold=0.5
fixed_contact_threshold=0
top_percent=50
min_raw_contacts=2
ccfind_terminal_fragment_size=500
ccfind_min_identity=94
ccfind_min_aligned_length=50
min_contact_strength=0
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
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

METAHICT_PATH=""
COMBINED=""
CONTACT_MATRIX=""
RAW_CONTACT_MATRIX=""
OUTDIR=""
THREADS=80
GENOMAD_DB=""
GENOMAD_SPLITS=8
GENOMAD_SENSITIVITY=4.2
GENOMAD_CLEANUP=true
GENOMAD_RESTART=false
GENOMAD_PRESET="default"
GENOMAD_MIN_SCORE=0.7
GENOMAD_MAX_FDR=0.1
GENOMAD_EXTRA_ARGS=""
ASSOCIATION_FILTER="zscore"
ZSCORE_THRESHOLD=0.5
FIXED_CONTACT_THRESHOLD=0
TOP_PERCENT=50
MIN_RAW_CONTACTS=2
CCFIND_TERMINAL_FRAGMENT_SIZE=500
CCFIND_MIN_IDENTITY=94
CCFIND_MIN_ALIGNED_LENGTH=50
MIN_CONTACT_STRENGTH=0
TMP_DIR="${METAHICT_TMP_ROOT:-${TMPDIR:-/tmp}}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--metahict-path|--project-path)
            METAHICT_PATH="$2"
            shift 2
            ;;
        --combined)
            COMBINED="$2"
            shift 2
            ;;
        --contact)
            CONTACT_MATRIX="$2"
            shift 2
            ;;
        --raw-contact)
            RAW_CONTACT_MATRIX="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --genomad_db|--genomad-db)
            GENOMAD_DB="$2"
            shift 2
            ;;
        --genomad-splits)
            GENOMAD_SPLITS="$2"
            shift 2
            ;;
        --genomad-sensitivity)
            GENOMAD_SENSITIVITY="$2"
            shift 2
            ;;
        --genomad-cleanup)
            GENOMAD_CLEANUP=true
            shift
            ;;
        --no-genomad-cleanup)
            GENOMAD_CLEANUP=false
            shift
            ;;
        --genomad-restart)
            GENOMAD_RESTART=true
            shift
            ;;
        --genomad-preset)
            GENOMAD_PRESET="$2"
            shift 2
            ;;
        --genomad-min-score)
            GENOMAD_MIN_SCORE="$2"
            shift 2
            ;;
        --genomad-max-fdr)
            GENOMAD_MAX_FDR="$2"
            shift 2
            ;;
        --genomad-extra-args)
            GENOMAD_EXTRA_ARGS="$2"
            shift 2
            ;;
        --association-filter)
            ASSOCIATION_FILTER="$2"
            shift 2
            ;;
        --zscore-threshold)
            ZSCORE_THRESHOLD="$2"
            shift 2
            ;;
        --fixed-contact-threshold)
            FIXED_CONTACT_THRESHOLD="$2"
            shift 2
            ;;
        --top-percent)
            TOP_PERCENT="$2"
            shift 2
            ;;
        --min-raw-contacts)
            MIN_RAW_CONTACTS="$2"
            shift 2
            ;;
        --ccfind-terminal-fragment-size)
            CCFIND_TERMINAL_FRAGMENT_SIZE="$2"
            shift 2
            ;;
        --ccfind-min-identity)
            CCFIND_MIN_IDENTITY="$2"
            shift 2
            ;;
        --ccfind-min-aligned-length)
            CCFIND_MIN_ALIGNED_LENGTH="$2"
            shift 2
            ;;
        --min-contact-strength)
            MIN_CONTACT_STRENGTH="$2"
            shift 2
            ;;
        --tmp-dir)
            TMP_DIR="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown parameter: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

[[ -n "$METAHICT_PATH" ]] || { echo "Error: --metahict-path is required" >&2; exit 1; }
[[ -f "$COMBINED" ]] || { echo "Error: combined FASTA not found: $COMBINED" >&2; exit 1; }
[[ -f "$CONTACT_MATRIX" ]] || { echo "Error: contact matrix not found: $CONTACT_MATRIX" >&2; exit 1; }
[[ -f "$RAW_CONTACT_MATRIX" ]] || { echo "Error: raw contact matrix not found: $RAW_CONTACT_MATRIX" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "Error: --outdir is required" >&2; exit 1; }
case "$GENOMAD_PRESET" in
    default|conservative|relaxed) ;;
    *) echo "Error: --genomad-preset must be default, conservative, or relaxed" >&2; exit 1 ;;
esac
case "$ASSOCIATION_FILTER" in
    zscore|fixed|percentage|raw-support-only) ;;
    *) echo "Error: --association-filter must be zscore, fixed, percentage, or raw-support-only" >&2; exit 1 ;;
esac
awk -v p="$TOP_PERCENT" 'BEGIN { exit !(p > 0 && p <= 100) }' || {
    echo "Error: --top-percent must be > 0 and <= 100" >&2
    exit 1
}
awk -v p="$MIN_RAW_CONTACTS" 'BEGIN { exit !(p >= 0) }' || {
    echo "Error: --min-raw-contacts must be >= 0" >&2
    exit 1
}

if [[ -z "$GENOMAD_DB" ]]; then
    GENOMAD_DB="$METAHICT_PATH/databases/genomad_db/genomad_db"
fi

mkdir -p "$OUTDIR/genomad_output" "$OUTDIR/mge_reports" "$TMP_DIR"
export TMPDIR="$TMP_DIR"

if command -v genomad >/dev/null 2>&1; then
    GENOMAD_BIN="$(command -v genomad)"
elif [[ -x "$METAHICT_PATH/conda_envs/genomad/bin/genomad" ]]; then
    GENOMAD_BIN="$METAHICT_PATH/conda_envs/genomad/bin/genomad"
else
    echo "Error: genomad executable not found." >&2
    exit 1
fi
export PATH="$(dirname "$GENOMAD_BIN"):$PATH"

if [[ -d "$METAHICT_PATH/conda_envs/ccfind_env/bin" ]]; then
    export PATH="$METAHICT_PATH/conda_envs/ccfind_env/bin:$PATH"
fi
if [[ -x "$METAHICT_PATH/external/bin/ccfind" ]]; then
    CCFIND_BIN="$METAHICT_PATH/external/bin/ccfind"
elif command -v ccfind >/dev/null 2>&1; then
    CCFIND_BIN="$(command -v ccfind)"
else
    echo "Error: ccfind executable not found." >&2
    exit 1
fi

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3)"
elif [[ -x "$METAHICT_PATH/conda_envs/metahict_env/bin/python" ]]; then
    PYTHON_BIN="$METAHICT_PATH/conda_envs/metahict_env/bin/python"
else
    echo "Error: Python executable not found." >&2
    exit 1
fi

[[ -d "$GENOMAD_DB" ]] || { echo "[ERROR] geNomad DB missing at $GENOMAD_DB" >&2; exit 1; }

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] Using geNomad executable: $GENOMAD_BIN"
echo "[INFO] Using ccfind executable: $CCFIND_BIN"

################################################################################
# Prepare a run-local combined FASTA with unique identifiers
################################################################################
echo "[INFO] ===== Preparing combined FASTA ====="
ORIGINAL_COMBINED="$COMBINED"
COMBINED_INPUT_DIR="$OUTDIR/input_combined"
mkdir -p "$COMBINED_INPUT_DIR"
COMBINED_UNIQUE="$COMBINED_INPUT_DIR/combined_contigs.unique.fa"
COMBINED_ID_MAP="$COMBINED_INPUT_DIR/combined_contigs.id_map.tsv"

export ORIGINAL_COMBINED COMBINED_UNIQUE COMBINED_ID_MAP
"$PYTHON_BIN" <<'EOF'
import os

source = os.environ["ORIGINAL_COMBINED"]
dest = os.environ["COMBINED_UNIQUE"]
map_file = os.environ["COMBINED_ID_MAP"]

seen = {}
used = set()
total = 0
renamed = 0

with open(source) as src, open(dest, "w") as out, open(map_file, "w") as mapping:
    mapping.write("original_id\tunique_id\toccurrence\n")
    for line in src:
        if not line.startswith(">"):
            out.write(line)
            continue
        total += 1
        header = line[1:].rstrip("\n")
        parts = header.split(maxsplit=1)
        original_id = parts[0]
        rest = f" {parts[1]}" if len(parts) > 1 else ""
        occurrence = seen.get(original_id, 0) + 1
        seen[original_id] = occurrence
        unique_id = original_id
        if occurrence > 1 or unique_id in used:
            renamed += 1
            suffix = occurrence
            unique_id = f"{original_id}__dup{suffix}"
            while unique_id in used:
                suffix += 1
                unique_id = f"{original_id}__dup{suffix}"
        used.add(unique_id)
        mapping.write(f"{original_id}\t{unique_id}\t{occurrence}\n")
        out.write(f">{unique_id}{rest}\n")

print(f"[INFO] Combined FASTA records: {total}")
print(f"[INFO] Duplicate FASTA IDs renamed: {renamed}")
print(f"[INFO] Combined FASTA ID map: {map_file}")
EOF
COMBINED="$COMBINED_UNIQUE"

################################################################################
# STEP 1: Run geNomad
################################################################################
echo "[INFO] ===== STEP 1: Run geNomad ====="
genomad_cmd=(
    "$GENOMAD_BIN" end-to-end
    --splits "$GENOMAD_SPLITS"
    --sensitivity "$GENOMAD_SENSITIVITY"
    --threads "$THREADS"
)
if [[ "$GENOMAD_CLEANUP" == true ]]; then
    genomad_cmd+=(--cleanup)
fi
if [[ "$GENOMAD_RESTART" == true ]]; then
    genomad_cmd+=(--restart)
fi
case "$GENOMAD_PRESET" in
    conservative) genomad_cmd+=(--conservative) ;;
    relaxed) genomad_cmd+=(--relaxed) ;;
    default)
        genomad_cmd+=(--min-score "$GENOMAD_MIN_SCORE" --max-fdr "$GENOMAD_MAX_FDR")
        ;;
esac
if [[ -n "$GENOMAD_EXTRA_ARGS" ]]; then
    read -r -a genomad_extra_args <<< "$GENOMAD_EXTRA_ARGS"
    genomad_cmd+=("${genomad_extra_args[@]}")
fi
genomad_cmd+=("$COMBINED" "$OUTDIR/genomad_output" "$GENOMAD_DB")

printf '[INFO] geNomad command:'
printf ' %q' "${genomad_cmd[@]}"
printf '\n'
"${genomad_cmd[@]}"
echo "[INFO] geNomad completed."

ASSEMBLY_BASE="$(basename "$COMBINED" .fa)"
ASSEMBLY_BASE="$(basename "$ASSEMBLY_BASE" .fasta)"
VIRUS_SUMMARY="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_virus_summary.tsv"
PLASMID_SUMMARY="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_plasmid_summary.tsv"
VIRAL_CONTIGS_RAW="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_virus.fna"
PLASMID_CONTIGS_RAW="$OUTDIR/genomad_output/${ASSEMBLY_BASE}_summary/${ASSEMBLY_BASE}_plasmid.fna"

################################################################################
# STEP 1.5: Remove proviruses
################################################################################
echo "[INFO] ===== STEP 1.5: Remove proviruses ====="
FILTERED_VIRAL="$OUTDIR/mge_reports/virus_no_provirus.fna"
mkdir -p "$OUTDIR/mge_reports"

export VIRUS_SUMMARY VIRAL_CONTIGS_RAW FILTERED_VIRAL
"$PYTHON_BIN" <<'EOF'
import os
import pandas as pd
from Bio import SeqIO

virus_summary = os.environ["VIRUS_SUMMARY"]
viral_fasta = os.environ["VIRAL_CONTIGS_RAW"]
filtered_out = os.environ["FILTERED_VIRAL"]

df = pd.read_csv(virus_summary, sep="\t")
keep_ids = set(df.loc[df["topology"] != "Provirus", "seq_name"])

with open(filtered_out, "w") as out_f:
    for rec in SeqIO.parse(viral_fasta, "fasta"):
        if rec.id in keep_ids:
            SeqIO.write(rec, out_f, "fasta")
EOF

VIRAL_QC="$FILTERED_VIRAL"

################################################################################
# STEP 2: Run ccfind on all combined contigs
################################################################################
echo "[INFO] ===== STEP 2: Run ccfind circularity detection ====="
CCFIND_OUTPUT_DIR="$OUTDIR/ccfind_output"
CCFIND_INPUT_DIR="$OUTDIR/ccfind_input"
CCFIND_INPUT="$CCFIND_INPUT_DIR/combined_contigs.ccfind_safe.fa"
CCFIND_ID_MAP="$CCFIND_INPUT_DIR/ccfind_id_map.tsv"
rm -rf "$CCFIND_OUTPUT_DIR"
mkdir -p "$CCFIND_INPUT_DIR"

export COMBINED CCFIND_INPUT CCFIND_ID_MAP
"$PYTHON_BIN" <<'EOF'
import os
from Bio import SeqIO

source = os.environ["COMBINED"]
dest = os.environ["CCFIND_INPUT"]
map_file = os.environ["CCFIND_ID_MAP"]

with open(dest, "w") as out, open(map_file, "w") as mapping:
    mapping.write("ccfind_id\tcontig_id\n")
    for idx, rec in enumerate(SeqIO.parse(source, "fasta"), 1):
        safe_id = f"ccfind_{idx:09d}"
        mapping.write(f"{safe_id}\t{rec.id}\n")
        rec.id = safe_id
        rec.name = safe_id
        rec.description = safe_id
        SeqIO.write(rec, out, "fasta")

print(f"[INFO] ccfind-safe FASTA written to: {dest}")
print(f"[INFO] ccfind ID map written to: {map_file}")
EOF

ccfind_cmd=(
    "$CCFIND_BIN"
    "$CCFIND_INPUT"
    "$CCFIND_OUTPUT_DIR"
    --terminal-fragment-size "$CCFIND_TERMINAL_FRAGMENT_SIZE"
    --min-percent-identity "$CCFIND_MIN_IDENTITY"
    --min-aligned-length "$CCFIND_MIN_ALIGNED_LENGTH"
)
if [[ "$THREADS" -gt 1 ]]; then
    ccfind_cmd+=(--ncpus "$THREADS")
fi
printf '[INFO] ccfind command:'
printf ' %q' "${ccfind_cmd[@]}"
printf '\n'
"${ccfind_cmd[@]}"
CCFIND_RESULT_DIR="$CCFIND_OUTPUT_DIR/result"
CCFIND_DETECTED_LIST_RAW="$CCFIND_RESULT_DIR/circ.detected.list"
CCFIND_DETECTED_LIST="$CCFIND_RESULT_DIR/circ.detected.mapped.list"
CCFIND_TOO_SHORT_LIST="$CCFIND_RESULT_DIR/too_short_seq.list"
touch "$CCFIND_DETECTED_LIST_RAW" "$CCFIND_TOO_SHORT_LIST"

export CCFIND_ID_MAP CCFIND_DETECTED_LIST_RAW CCFIND_DETECTED_LIST
"$PYTHON_BIN" <<'EOF'
import os

map_file = os.environ["CCFIND_ID_MAP"]
raw_file = os.environ["CCFIND_DETECTED_LIST_RAW"]
mapped_file = os.environ["CCFIND_DETECTED_LIST"]

safe_to_contig = {}
with open(map_file) as handle:
    next(handle, None)
    for line in handle:
        safe_id, contig_id = line.rstrip("\n").split("\t", 1)
        safe_to_contig[safe_id] = contig_id

with open(raw_file) as inp, open(mapped_file, "w") as out:
    for line in inp:
        stripped = line.rstrip("\n")
        if not stripped:
            continue
        fields = stripped.split("\t")
        fields[0] = safe_to_contig.get(fields[0], fields[0])
        out.write("\t".join(fields) + "\n")

print(f"[INFO] Mapped ccfind circular IDs written to: {mapped_file}")
EOF
echo "[INFO] ccfind completed."

################################################################################
# Reports and linkages
################################################################################
echo "[INFO] ===== Reporting ====="

export OUTDIR VIRUS_SUMMARY PLASMID_SUMMARY FILTERED_VIRAL PLASMID_CONTIGS_RAW
export VIRAL_QC CONTACT_MATRIX RAW_CONTACT_MATRIX COMBINED
export ASSOCIATION_FILTER ZSCORE_THRESHOLD FIXED_CONTACT_THRESHOLD TOP_PERCENT MIN_RAW_CONTACTS
export CCFIND_RESULT_DIR CCFIND_DETECTED_LIST CCFIND_TOO_SHORT_LIST CCFIND_TERMINAL_FRAGMENT_SIZE CCFIND_MIN_IDENTITY CCFIND_MIN_ALIGNED_LENGTH MIN_CONTACT_STRENGTH

"$PYTHON_BIN" <<'EOF'
import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
try:
    from scipy.sparse import load_npz as scipy_load_npz
except ModuleNotFoundError:
    scipy_load_npz = None

class SimpleCOO:
    def __init__(self, row, col, data, shape):
        self.row = row
        self.col = col
        self.data = data
        self.shape = tuple(int(x) for x in shape)

def load_npz_coo(path):
    if scipy_load_npz is not None:
        return scipy_load_npz(path).tocoo()

    archive = np.load(path, allow_pickle=False)
    fmt = archive["format"]
    if hasattr(fmt, "item"):
        fmt = fmt.item()
    if isinstance(fmt, bytes):
        fmt = fmt.decode()
    fmt = str(fmt)
    shape = tuple(int(x) for x in archive["shape"])
    data = archive["data"]

    if fmt == "coo":
        return SimpleCOO(archive["row"], archive["col"], data, shape)
    if fmt == "csr":
        indptr = archive["indptr"]
        indices = archive["indices"]
        row = np.repeat(np.arange(shape[0]), np.diff(indptr))
        return SimpleCOO(row, indices, data, shape)
    if fmt == "csc":
        indptr = archive["indptr"]
        indices = archive["indices"]
        col = np.repeat(np.arange(shape[1]), np.diff(indptr))
        return SimpleCOO(indices, col, data, shape)
    raise ValueError(f"Unsupported sparse matrix format in {path}: {fmt}")

outdir = os.environ["OUTDIR"]
virus_summary = os.environ["VIRUS_SUMMARY"]
plasmid_summary = os.environ["PLASMID_SUMMARY"]
filtered_viral = os.environ["FILTERED_VIRAL"]
plasmid_raw = os.environ["PLASMID_CONTIGS_RAW"]
viral_qc = os.environ["VIRAL_QC"]
contact_file = os.environ["CONTACT_MATRIX"]
raw_contact_file = os.environ["RAW_CONTACT_MATRIX"]
combined = os.environ["COMBINED"]
association_filter = os.environ.get("ASSOCIATION_FILTER", "zscore")
zscore_threshold = float(os.environ.get("ZSCORE_THRESHOLD", "0.5"))
fixed_contact_threshold = float(os.environ.get("FIXED_CONTACT_THRESHOLD", "0"))
top_percent = float(os.environ.get("TOP_PERCENT", "50"))
min_raw_contacts = float(os.environ.get("MIN_RAW_CONTACTS", "2"))
ccfind_result_dir = os.environ["CCFIND_RESULT_DIR"]
ccfind_detected_list = os.environ["CCFIND_DETECTED_LIST"]
ccfind_too_short_list = os.environ["CCFIND_TOO_SHORT_LIST"]
ccfind_terminal_fragment_size = int(os.environ.get("CCFIND_TERMINAL_FRAGMENT_SIZE", "500"))
ccfind_min_identity = int(os.environ.get("CCFIND_MIN_IDENTITY", "94"))
ccfind_min_aligned_length = int(os.environ.get("CCFIND_MIN_ALIGNED_LENGTH", "50"))
min_contact_strength = float(os.environ.get("MIN_CONTACT_STRENGTH", "0"))

summary_file = f"{outdir}/mge_reports/MGE_summary.txt"

def count_bin_unmapped(fasta):
    bins = sum(1 for rec in SeqIO.parse(fasta, "fasta") if rec.id.startswith("bin"))
    total = sum(1 for _ in SeqIO.parse(fasta, "fasta"))
    unmapped = total - bins
    return total, bins, unmapped

def fasta_lengths(fasta):
    if not os.path.exists(fasta):
        return {}
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}

def fasta_sequences(fasta):
    if not os.path.exists(fasta):
        return {}
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta, "fasta")}

def source_from_contig(contig_id):
    return "bin" if contig_id.startswith("bin") else "unmapped"

def read_ccfind_detected(path):
    detected = set()
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return detected
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            detected.add(line.split("\t", 1)[0])
    return detected

def count_nonempty_lines(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0
    count = 0
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count

def read_summary(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()

def summary_id_column(df):
    for col in ("seq_name", "contig_id", "id"):
        if col in df.columns:
            return col
    return df.columns[0] if len(df.columns) else None

def topology_map(df):
    id_col = summary_id_column(df)
    if id_col is None or "topology" not in df.columns:
        return {}
    return dict(zip(df[id_col].astype(str), df["topology"].fillna("").astype(str)))

def summary_ids(df):
    id_col = summary_id_column(df)
    if id_col is None:
        return set()
    return set(df[id_col].astype(str))

def genomad_circular_evidence(topology):
    topology_text = str(topology).lower()
    if any(term in topology_text for term in ("circular", "dtr", "itr")):
        return "yes"
    return "no"

def host_mag_from_contig(contig_id):
    if not contig_id.startswith("bin"):
        return contig_id
    for sep in ("|", ":", ";"):
        if sep in contig_id:
            return contig_id.split(sep, 1)[0]
    match = re.match(r"^(bin[^_]+)", contig_id)
    return match.group(1) if match else contig_id

def extract_taxonomy_map(df):
    if len(df.columns) == 0:
        return {}
    id_col = "seq_name" if "seq_name" in df.columns else df.columns[0]
    tax_cols = [
        col for col in df.columns
        if any(token in col.lower() for token in ("tax", "lineage", "virus_name", "classification"))
    ]
    if not tax_cols:
        return {}
    tax_df = df[[id_col] + tax_cols].copy()
    tax_df["mge_taxonomy"] = tax_df[tax_cols].fillna("").astype(str).agg(";".join, axis=1)
    return dict(zip(tax_df[id_col], tax_df["mge_taxonomy"]))

def contact_contigs_for_matrix(matrix_file, combined_fasta):
    contact_dir = os.path.dirname(os.path.abspath(matrix_file))
    contig_info = os.path.join(contact_dir, "contig_info.csv")
    if os.path.exists(contig_info):
        info = pd.read_csv(contig_info)
        if "name" not in info.columns:
            raise ValueError(f"Contact contig info file lacks a 'name' column: {contig_info}")
        return info["name"].astype(str).tolist(), contig_info
    contigs = [line[1:].split()[0] for line in open(combined_fasta) if line.startswith(">")]
    return contigs, combined_fasta

with open(summary_file, "w") as fh:
    fh.write("=== STEP 1.5 Report (After Provirus Removal) ===\n")
    df = pd.read_csv(virus_summary, sep="\t")
    provirus_count = (df["topology"] == "Provirus").sum()
    virus_total, virus_bins, virus_unmapped = count_bin_unmapped(filtered_viral)
    plasmid_total, plasmid_bins, plasmid_unmapped = count_bin_unmapped(plasmid_raw)
    fh.write(f"Virus contigs (after provirus removal): {virus_total}\n")
    fh.write(f"  from bins: {virus_bins}\n")
    fh.write(f"  from unmapped: {virus_unmapped}\n")
    fh.write(f"Provirus removed: {provirus_count}\n")
    fh.write(f"Plasmid contigs (raw): {plasmid_total}\n")
    fh.write(f"  from bins: {plasmid_bins}\n")
    fh.write(f"  from unmapped: {plasmid_unmapped}\n\n")

    fh.write("=== STEP 2 Report (Viral Candidates for Host Attribution) ===\n")
    ids_virus = set(fasta_lengths(filtered_viral))
    fh.write(f"Virus total (for host attribution): {len(ids_virus)}\n\n")

    fh.write("=== STEP 2.5 Report (Sequence Topology) ===\n")
    plasmid_df = read_summary(plasmid_summary)
    virus_topology = topology_map(df)
    plasmid_topology = topology_map(plasmid_df)
    viral_lengths = fasta_lengths(filtered_viral)
    plasmid_lengths = fasta_lengths(plasmid_raw)
    ccfind_detected = read_ccfind_detected(ccfind_detected_list)
    ccfind_too_short_count = count_nonempty_lines(ccfind_too_short_list)
    plasmid_ids = set(plasmid_lengths)
    combined_sequences = fasta_sequences(combined)
    mge_ids = ids_virus | plasmid_ids

    sequence_topology_rows = []

    for contig_id in sorted(combined_sequences):
        is_mge = contig_id in mge_ids
        ccfind = 1 if contig_id in ccfind_detected else 0
        genomad = "NA"
        if contig_id in ids_virus:
            topology = virus_topology.get(contig_id, "")
            genomad = 1 if genomad_circular_evidence(topology) == "yes" else 0
        elif contig_id in plasmid_ids:
            topology = plasmid_topology.get(contig_id, "")
            genomad = 1 if genomad_circular_evidence(topology) == "yes" else 0

        sequence_topology_rows.append({
            "contig_id": contig_id,
            "mge": 1 if is_mge else 0,
            "ccfind": ccfind,
            "genomad": genomad,
        })

    sequence_topology_file = f"{outdir}/sequence_topology.tsv"
    sequence_topology_columns = ["contig_id", "mge", "ccfind", "genomad"]
    sequence_topology_df = pd.DataFrame(sequence_topology_rows).reindex(columns=sequence_topology_columns)
    circular_viruses = int(((sequence_topology_df["contig_id"].isin(ids_virus)) & ((sequence_topology_df["ccfind"] == 1) | (sequence_topology_df["genomad"] == 1))).sum())
    circular_plasmids = int(((sequence_topology_df["contig_id"].isin(plasmid_ids)) & ((sequence_topology_df["ccfind"] == 1) | (sequence_topology_df["genomad"] == 1))).sum())
    genomad_circular_mge_count = int((sequence_topology_df["genomad"] == 1).sum())
    ccfind_circular_mge_count = int(((sequence_topology_df["mge"] == 1) & (sequence_topology_df["ccfind"] == 1)).sum())
    ccfind_circular_non_mge_count = int(((sequence_topology_df["mge"] == 0) & (sequence_topology_df["ccfind"] == 1)).sum())
    pd.DataFrame(sequence_topology_rows).reindex(columns=sequence_topology_columns).to_csv(
        sequence_topology_file, sep="\t", index=False
    )
    fh.write(
        "ccfind parameters: "
        f"terminal_fragment_size={ccfind_terminal_fragment_size}, "
        f"min_percent_identity={ccfind_min_identity}, "
        f"min_aligned_length={ccfind_min_aligned_length}\n"
    )
    fh.write(f"ccfind too-short contigs excluded: {ccfind_too_short_count}\n")
    fh.write(f"Combined contigs evaluated: {len(sequence_topology_rows)}\n")
    fh.write(f"Circular viral candidates by geNomad or ccfind: {circular_viruses}\n")
    fh.write(f"Circular plasmid candidates by geNomad or ccfind: {circular_plasmids}\n")
    fh.write(f"Circular MGE candidates by geNomad topology: {genomad_circular_mge_count}\n")
    fh.write(f"Circular MGE candidates by ccfind: {ccfind_circular_mge_count}\n")
    fh.write(f"ccfind circular non-MGE contigs: {ccfind_circular_non_mge_count}\n")
    fh.write("Sequence topology report saved to: sequence_topology.tsv\n\n")

    fh.write("=== STEP 3 Report (Candidate MGE-MAG Contacts) ===\n")
    normalized_mat = load_npz_coo(contact_file)
    raw_mat = load_npz_coo(raw_contact_file)
    contigs, contig_order_source = contact_contigs_for_matrix(contact_file, combined)
    raw_contigs, raw_contig_order_source = contact_contigs_for_matrix(raw_contact_file, combined)
    if raw_contigs != contigs:
        raise ValueError("Raw and normalized contact matrices use different contig orders")
    if normalized_mat.shape[0] != len(contigs) or normalized_mat.shape[1] != len(contigs):
        raise ValueError(
            f"Contact matrix shape {normalized_mat.shape} does not match contact contig count {len(contigs)}"
        )
    if raw_mat.shape[0] != len(raw_contigs) or raw_mat.shape[1] != len(raw_contigs):
        raise ValueError(
            f"Raw contact matrix shape {raw_mat.shape} does not match contact contig count {len(raw_contigs)}"
        )
    viral_contigs = set([line[1:].split()[0] for line in open(viral_qc) if line.startswith(">")])
    plasmid_contigs = set([line[1:].split()[0] for line in open(plasmid_raw) if line.startswith(">")])
    mge_type_map = {contig: "virus" for contig in viral_contigs}
    mge_type_map.update({contig: "plasmid" for contig in plasmid_contigs if contig not in mge_type_map})
    mge_contigs = set(mge_type_map)
    host_contigs = {c for c in contigs if (c.startswith("bin") and c not in mge_contigs)}
    def count_contacts(mge_set):
        count = 0
        for r, c, v in zip(normalized_mat.row, normalized_mat.col, normalized_mat.data):
            if r >= c:
                continue
            if float(v) <= min_contact_strength:
                continue
            a, b = contigs[r], contigs[c]
            if (a in mge_set and b in host_contigs) or (b in mge_set and a in host_contigs):
                count += 1
        return count

    fh.write("Association score source: normalized\n")
    fh.write(f"Association normalized contact matrix: {contact_file}\n")
    fh.write(f"Association raw contact matrix: {raw_contact_file}\n")
    fh.write(f"Association contact contig order: {contig_order_source}\n")
    fh.write(f"Association raw contact contig order: {raw_contig_order_source}\n")
    fh.write(f"Association filter: {association_filter}\n")
    fh.write(f"Minimum normalized contact strength: {min_contact_strength}\n")
    fh.write(f"Minimum raw Hi-C support: {min_raw_contacts}\n")
    fh.write(f"Viral contigs: {len(viral_contigs)}\n")
    fh.write(f"Plasmid contigs: {len(plasmid_contigs)}\n")
    fh.write(f"Viral candidate MGE-MAG contacts before Z-score filtering: {count_contacts(viral_contigs)}\n")
    fh.write(f"Plasmid candidate MGE-MAG contacts before Z-score filtering: {count_contacts(plasmid_contigs)}\n\n")

    fh.write("=== STEP 4 Report (Candidate MGE-MAG Association Table) ===\n")
    taxonomy_map = extract_taxonomy_map(df)
    taxonomy_map.update(extract_taxonomy_map(plasmid_df))
    normalized_associations = {}
    for r, c, v in zip(normalized_mat.row, normalized_mat.col, normalized_mat.data):
        if r >= c:
            continue
        if float(v) <= min_contact_strength:
            continue
        contig_a = contigs[r]
        contig_b = contigs[c]
        if contig_a in mge_contigs and contig_b in host_contigs:
            mge_contig = contig_a
            host_contig = contig_b
        elif contig_b in mge_contigs and contig_a in host_contigs:
            mge_contig = contig_b
            host_contig = contig_a
        else:
            continue

        host_mag = host_mag_from_contig(host_contig)
        key = (mge_contig, host_mag)
        if key not in normalized_associations:
            normalized_associations[key] = {
                "mge_contig": mge_contig,
                "mge_type": mge_type_map.get(mge_contig, ""),
                "host_mag": host_mag,
                "_host_contigs": set(),
                "normalized_contact_score": 0.0,
                "max_contig_normalized_contact_score": 0.0,
                "mge_taxonomy": taxonomy_map.get(mge_contig, "")
            }
        normalized_associations[key]["_host_contigs"].add(host_contig)
        normalized_associations[key]["normalized_contact_score"] += float(v)
        normalized_associations[key]["max_contig_normalized_contact_score"] = max(
            normalized_associations[key]["max_contig_normalized_contact_score"], float(v)
        )

    raw_support = {}
    for r, c, v in zip(raw_mat.row, raw_mat.col, raw_mat.data):
        if r >= c:
            continue
        if float(v) <= 0:
            continue
        contig_a = raw_contigs[r]
        contig_b = raw_contigs[c]
        if contig_a in mge_contigs and contig_b in host_contigs:
            mge_contig = contig_a
            host_contig = contig_b
        elif contig_b in mge_contigs and contig_a in host_contigs:
            mge_contig = contig_b
            host_contig = contig_a
        else:
            continue
        host_mag = host_mag_from_contig(host_contig)
        key = (mge_contig, host_mag)
        raw_support[key] = raw_support.get(key, 0.0) + float(v)

    for key, association in normalized_associations.items():
        association["host_contig_count"] = len(association["_host_contigs"])
        del association["_host_contigs"]
        association["raw_contact_support"] = raw_support.get(key, 0.0)

    association_columns = [
        "mge_contig", "mge_type", "host_mag", "raw_contact_support",
        "normalized_contact_score", "max_contig_normalized_contact_score",
        "mge_taxonomy", "host_contig_count", "z_score",
        "passes_raw_contact_filter",
        "passes_association_filter", "passes_zscore_threshold",
        "association_score_source", "association_filter"
    ]
    filtered_file = f"{outdir}/candidate_mge_mag_associations_{association_filter}_filtered.tsv"
    associations_df = pd.DataFrame(normalized_associations.values())
    if len(associations_df) > 0:
        associations_df["z_score"] = 0.0
        for _, idx in associations_df.groupby("mge_type").groups.items():
            scores = associations_df.loc[idx, "normalized_contact_score"]
            mean = scores.mean()
            std = scores.std(ddof=0)
            if std != 0 and not pd.isna(std):
                associations_df.loc[idx, "z_score"] = (scores - mean) / std
        associations_df["passes_raw_contact_filter"] = associations_df["raw_contact_support"] >= min_raw_contacts
        if association_filter == "zscore":
            associations_df["passes_association_filter"] = (
                associations_df["passes_raw_contact_filter"] &
                (associations_df["z_score"] > zscore_threshold)
            )
        elif association_filter == "fixed":
            associations_df["passes_association_filter"] = (
                associations_df["passes_raw_contact_filter"] &
                (associations_df["normalized_contact_score"] >= fixed_contact_threshold)
            )
        elif association_filter == "percentage":
            n_keep = max(1, int((len(associations_df) * top_percent + 99.999999) // 100))
            associations_df["_rank"] = associations_df["normalized_contact_score"].rank(method="first", ascending=False)
            associations_df["passes_association_filter"] = (
                associations_df["passes_raw_contact_filter"] &
                (associations_df["_rank"] <= n_keep)
            )
            associations_df = associations_df.drop(columns=["_rank"])
        elif association_filter == "raw-support-only":
            associations_df["passes_association_filter"] = associations_df["passes_raw_contact_filter"]
        else:
            raise ValueError(f"Unsupported association filter: {association_filter}")
        associations_df["passes_zscore_threshold"] = associations_df["z_score"] > zscore_threshold
        associations_df["association_score_source"] = "normalized"
        associations_df["association_filter"] = association_filter
        associations_df = associations_df.sort_values(
            ["passes_association_filter", "z_score", "normalized_contact_score"],
            ascending=[False, False, False]
        )
        associations_df = associations_df.reindex(columns=association_columns)

        filtered_df = associations_df[associations_df["passes_association_filter"]]
        filtered_df.to_csv(filtered_file, sep="\t", index=False)

        fh.write("Association score source: normalized\n")
        fh.write(f"Association filter: {association_filter}\n")
        fh.write(f"Z-score threshold: {zscore_threshold}\n")
        fh.write(f"Fixed contact threshold: {fixed_contact_threshold}\n")
        fh.write(f"Top percent: {top_percent}\n")
        fh.write(f"Minimum raw Hi-C support: {min_raw_contacts}\n")
        fh.write(f"Scored candidate MGE-MAG associations before filtering: {len(associations_df)}\n")
        fh.write(f"Filtered candidate MGE-MAG associations: {len(filtered_df)}\n")
        fh.write(f"Unique MGE contigs after filtering: {filtered_df['mge_contig'].nunique()}\n")
        for mge_type, sub_df in filtered_df.groupby("mge_type"):
            fh.write(f"Unique {mge_type} MGE contigs after filtering: {sub_df['mge_contig'].nunique()}\n")
        fh.write(f"Unique host MAGs after filtering: {filtered_df['host_mag'].nunique()}\n")
        fh.write(f"Filtered candidate associations table saved to: {os.path.basename(filtered_file)}\n")
    else:
        pd.DataFrame(columns=association_columns).to_csv(filtered_file, sep="\t", index=False)
        fh.write("No candidate MGE-MAG associations found.\n")
        fh.write(f"Filtered candidate associations table saved to: {os.path.basename(filtered_file)}\n")

EOF

echo "[INFO] ===== MGE analysis finished ====="
echo "[INFO] Summary: $OUTDIR/mge_reports/MGE_summary.txt"
echo "[INFO] Candidate associations: $OUTDIR/candidate_mge_mag_associations_${ASSOCIATION_FILTER}_filtered.tsv"
