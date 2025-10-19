#!/bin/bash

# Define default download directory
DB_DIR="$../databases"
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

# 2. Download CheckM2 Database
checkm2 database --download --path "$DB_DIR"

# Cleanup downloaded files
rm -f "$CHECKM2_DB" 

echo "All databases downloaded and extracted to: $DB_DIR"