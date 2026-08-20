#!/bin/bash
set -euo pipefail

# genomad_db.sh
# Download and set up the geNomad database for MetaHit.
# By default, installs into "databases/genomad_db" under the MetaHit root.

# Determine database directory
if [ -n "${1-}" ]; then
    DB_DIR="$1"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DB_DIR="$(cd "$SCRIPT_DIR/../../" && pwd)/databases/genomad_db"
fi

mkdir -p "$DB_DIR"
echo "[INFO] Database directory set to: $DB_DIR"

# Function to print messages
echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

# This script should be run from an active environment that provides genomad.
# It does not activate conda internally because some activation hooks are not
# compatible with strict shell mode.
if ! command -v genomad &> /dev/null; then
    echo "[ERROR] geNomad executable not found."
    echo "[ERROR] Please run: conda activate genomad"
    exit 1
fi

# Download geNomad database
echo_info "Downloading geNomad database to $DB_DIR..."
genomad download-database "$DB_DIR"

# Verify database structure
if [ -f "$DB_DIR/genomad_db/version.txt" ]; then
    echo_info "geNomad database successfully downloaded."
else
    echo "[ERROR] geNomad database download failed or incomplete."
    exit 1
fi

# Output success message
echo_info "geNomad database setup completed successfully."
echo_info "Database installed at: $DB_DIR"
echo_info "You can later set METAHICT_PATH/databases/genomad_db to point here if needed."
