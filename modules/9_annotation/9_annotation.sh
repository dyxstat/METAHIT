#!/usr/bin/env bash

# Print free memory
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
echo "Running: gtdbtk classify_wf"

# Display usage if insufficient arguments
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 -p <METAHIT_PATH> --bin <BIN_DIR> --outdir <OUTPUT_DIR> -t <THREADS> [--gtdbtk_db <GTDB_DB_PATH>] [additional gtdbtk options]"
    exit 1
fi

# Defaults
THREADS=1
GTDBTK_DB=""

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p)
            PATH_DIR=$2
            shift 2
            ;;
        --bin)
            BIN_DIR=$2
            shift 2
            ;;
        --outdir)
            OUT_DIR=$2
            shift 2
            ;;
        -t)
            THREADS=$2
            shift 2
            ;;
        --gtdbtk_db)
            GTDBTK_DB=$2
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

# Validate required arguments
if [ -z "$PATH_DIR" ] || [ -z "$BIN_DIR" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Missing required argument(s)."
    exit 1
fi

# Activate the GTDB-Tk conda environment
echo "[INFO] Activating GTDB-Tk environment 'gtdbtk-2.4.0'..."
eval "$(conda shell.bash hook)"
conda activate gtdbtk-2.4.0

if [ $? -ne 0 ]; then
    echo "Error: Could not activate Conda environment 'gtdbtk-2.4.0'."
    exit 1
fi

# Set GTDB-Tk database path
if [ -n "$GTDBTK_DB" ]; then
    export GTDBTK_DATA_PATH="$GTDBTK_DB"
    echo "[INFO] Using custom GTDB-Tk database: $GTDBTK_DB"
else
    export GTDBTK_DATA_PATH="${PATH_DIR}/databases/release220"
    echo "[INFO] Using default GTDB-Tk database: ${PATH_DIR}/databases/release220"
fi

# Run GTDB-Tk classification
echo "[INFO] Running GTDB-Tk classify_wf with $THREADS threads..."
gtdbtk classify_wf \
    --genome_dir "$BIN_DIR" \
    --out_dir "$OUT_DIR" \
    --extension fa \
    --cpus "$THREADS" \
    --pplacer_cpus "$THREADS" \
    --skip_ani_screen "$@"

if [ $? -ne 0 ]; then
    echo "Error: Annotation (GTDB-Tk) step failed."
    exit 1
fi

conda deactivate
conda activate metahit_env

echo "Annotation (GTDB-Tk) step completed successfully."
