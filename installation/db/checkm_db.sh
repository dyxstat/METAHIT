#!/bin/bash

# checkm_db.sh
# Download and extract the CheckM database to a user-specified or default directory.

# If user provides a directory, use it; otherwise default to two levels up (METAHICT/databases)
if [ -n "$1" ]; then
    DB_DIR="$1"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DB_DIR="$(cd "$SCRIPT_DIR/../../" && pwd)/databases"
fi

mkdir -p "$DB_DIR"

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

# 1. Download CheckM Database
CHECKM_URL="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
CHECKM_DB="$DB_DIR/checkm_data.tar.gz"

download_file "$CHECKM_URL" "$CHECKM_DB"
extract_tar_gz "$CHECKM_DB" "$DB_DIR"

# Set CheckM database root
checkm data setRoot "$DB_DIR"

# Cleanup downloaded tarball
rm -f "$CHECKM_DB"

echo "[INFO] CheckM database installed to: $DB_DIR"
