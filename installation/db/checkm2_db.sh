#!/bin/bash

# checkm2_db.sh
# Download and set up the CheckM2 database in a user-specified or default directory.

# Determine database directory
if [ -n "$1" ]; then
    DB_DIR="$1"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DB_DIR="$(cd "$SCRIPT_DIR/../../" && pwd)/databases"
fi

mkdir -p "$DB_DIR"

echo "[INFO] Database directory set to: $DB_DIR"

# Run CheckM2 database download
echo "[INFO] Downloading CheckM2 database..."
checkm2 database --download --path "$DB_DIR"

echo "[INFO] CheckM2 database successfully installed at: $DB_DIR"
