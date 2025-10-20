#!/bin/bash
set -euo pipefail

# gtdbtk_db.sh
# Download and set up the GTDB-Tk database for MetaHit.

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

# Download GTDB-Tk Database
GTDBTK_URL="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
GTDBTK_DB="$DB_DIR/gtdbtk_data.tar.gz"

download_file "$GTDBTK_URL" "$GTDBTK_DB"
extract_tar_gz "$GTDBTK_DB" "$DB_DIR"

# Configure environment variable for GTDB-Tk
echo "[INFO] Setting GTDBTK_DATA_PATH..."
GTDBTK_PATH_ABS="$(realpath "$DB_DIR/release220")"
conda env config vars set GTDBTK_DATA_PATH="$GTDBTK_PATH_ABS"

# Cleanup compressed archive
echo "[INFO] Cleaning up compressed archive..."
rm -f "$GTDBTK_DB"

echo "[INFO] GTDB-Tk database successfully installed at: $DB_DIR"
