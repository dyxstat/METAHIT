#!/bin/bash
set -euo pipefail

# checkv_db.sh
# Download and set up the CheckV database for METAHICT.

# Determine database directory
if [ -n "${1-}" ]; then
    DB_DIR="$1"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DB_DIR="$(cd "$SCRIPT_DIR/../../" && pwd)/databases"
fi

mkdir -p "$DB_DIR"

echo "[INFO] Database directory set to: $DB_DIR"

# Function to download a file with progress
download_file() {
    local url=$1
    local dest=$2
    echo "[INFO] Downloading: $url"
    wget -c "$url" -O "$dest"
}

# Function to extract tar.gz files
extract_tar_gz() {
    local file=$1
    local dest_dir=$2
    echo "[INFO] Extracting: $file"
    tar -xvzf "$file" -C "$dest_dir"
}

# Download CheckV Database
CHECKV_URL="https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
CHECKV_DB="$DB_DIR/checkv-db-v1.5.tar.gz"

download_file "$CHECKV_URL" "$CHECKV_DB"
extract_tar_gz "$CHECKV_DB" "$DB_DIR"

if command -v diamond >/dev/null 2>&1; then
    if [ -f "$DB_DIR/checkv-db-v1.5/genome_db/checkv_reps.faa" ] && [ ! -f "$DB_DIR/checkv-db-v1.5/genome_db/checkv_reps.dmnd" ]; then
        echo "[INFO] Building CheckV DIAMOND database..."
        diamond makedb \
            --in "$DB_DIR/checkv-db-v1.5/genome_db/checkv_reps.faa" \
            --db "$DB_DIR/checkv-db-v1.5/genome_db/checkv_reps"
    fi
else
    echo "[WARNING] diamond was not found in PATH; CheckV may build or require checkv_reps.dmnd at runtime."
fi

# Cleanup compressed archive
echo "[INFO] Cleaning up compressed archive..."
rm -f "$CHECKV_DB"

echo "[INFO] CheckV database successfully installed at: $DB_DIR"
