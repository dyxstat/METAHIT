#!/bin/bash

set -euo pipefail

# Define default download directory
DB_DIR="$(pwd)/databases"
mkdir -p "$DB_DIR"

# Function to download a file with progress
download_file() {
    local url=$1
    local dest=$2
    echo "Downloading: $url"
    wget -c "$url" -O "$dest"
}

# Function to extract tar.gz files
extract_tar_gz() {
    local file=$1
    local dest_dir=$2
    echo "Extracting: $file"
    tar -xvzf "$file" -C "$dest_dir"
}

# 3. Download GTDB-Tk Database
GTDBTK_URL="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
GTDBTK_DB="$DB_DIR/gtdbtk_data.tar.gz"

download_file "$GTDBTK_URL" "$GTDBTK_DB"
extract_tar_gz "$GTDBTK_DB" "$DB_DIR"

echo "Activating environment and setting GTDBTK_DATA_PATH..."
GTDBTK_PATH_ABS="$(realpath "$DB_DIR/release220")"
conda env config vars set GTDBTK_DATA_PATH="$GTDBTK_PATH_ABS"

echo "Cleaning up compressed archive..."
rm -f "$GTDBTK_DB"

echo "All databases downloaded and extracted to: $DB_DIR"