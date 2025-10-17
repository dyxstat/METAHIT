#!/bin/bash
set -euo pipefail

# checkv_db.sh
# Download and set up the CheckV database for MetaHit.

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

# 4. Download CheckV Database
CHECKV_URL="https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
CHECKV_DB="$DB_DIR/checkv-db-v1.5.tar.gz"

download_file "$CHECKV_URL" "$CHECKV_DB"
extract_tar_gz "$CHECKV_DB" "$DB_DIR"

# Cleanup compressed archive
echo "Cleaning up compressed archive..."
rm -f "$CHECKV_DB"

echo "All databases downloaded and extracted to: $DB_DIR"
