#!/bin/bash

# setup.sh
# This script sets up the necessary dependencies for the MetaHit pipeline.
# into the "external" directory within the repository.

# Exit immediately if a command exits with a non-zero status
set -e

# Function to print informational messages
function echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

# Function to print error messages
function echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

# Define the external and bin directory paths
EXTERNAL_DIR="$(pwd)/external"
BIN_DIR="${EXTERNAL_DIR}/bin"

# Create the external and bin directories if they don't exist
if [ ! -d "$EXTERNAL_DIR" ]; then
    echo_info "Creating 'external' directory."
    mkdir -p "$EXTERNAL_DIR"
else
    echo_info "'external' directory already exists."
fi

if [ ! -d "$BIN_DIR" ]; then
    echo_info "Creating 'bin' directory inside 'external'."
    mkdir -p "$BIN_DIR"
else
    echo_info "'bin' directory already exists inside 'external'."
fi

# Install dependencies via Conda
echo_info "Installing dependencies using Conda..."
conda install -y -c bioconda wget unzip openjdk perl git



function install_bbtools() {
    BBTOOLS_VERSION="39.10"  # Latest version as per user request
    BBTOOLS_TARBALL="Bbmap_${BBTOOLS_VERSION}.tar.gz"
    BBTOOLS_URL="https://sourceforge.net/projects/bbmap/files/latest/download"
    if [ ! -f "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" ]; then
        echo_info "Downloading BBTools ${BBTOOLS_VERSION}..."
        wget -O "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" "${BBTOOLS_URL}"
    else
        echo_info "BBTools tarball already exists, skipping download."
    fi

    if [ ! -f "${BIN_DIR}/bbmap.sh" ]; then
        echo_info "Extracting BBTools..."
        tar -xzf "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" -C "$EXTERNAL_DIR" || { echo_error "Failed to extract BBTools tarball."; exit 1; }

        echo_info "Creating symbolic links for BBTools binaries in 'external/bin'..."
        cd "${EXTERNAL_DIR}/bbmap"
        for bin in *.sh; do
            ln -sf "${EXTERNAL_DIR}/bbmap/${bin}" "${BIN_DIR}/${bin}"
        done

        cd ../..

        echo_info "BBTools installed successfully."
    else
        echo_info "BBTools binaries already exist in 'external/bin', skipping extraction."
    fi
}

# Install all dependencies
install_bbtools 
# Verify installations
echo_info "Verifying installations..."



# Check if the gtdbtk-2.4.0 environment exists
if ! conda info --envs | grep -q "gtdbtk-2.4.0"; then
    echo "[INFO] Creating GTDB-Tk environment 'gtdbtk-2.4.0'..."
    conda create -n gtdbtk-2.4.0 -c bioconda -c conda-forge gtdbtk=2.4.0 -y
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create GTDB-Tk environment."
        exit 1
    fi
fi


conda env create -f env.yaml
conda env create -f checkm2.yaml
conda create -n genomad -c conda-forge -c bioconda genomad

# Create CheckV environment
echo_info "Creating CheckV environment 'checkv_env'..."
if ! conda info --envs | grep -q "checkv_env"; then
    conda create -n checkv_env -y -c bioconda -c conda-forge \
        checkv \
        python=3.10 \
        biopython \
        pandas \
        numpy \
        h5py \
        scipy \
        prodigal \
        hmmer \
        diamond \
        tqdm
    if [ $? -ne 0 ]; then
        echo_error "Failed to create CheckV environment."
        exit 1
    fi
else
    echo_info "Environment 'checkv_env' already exists. Skipping creation."
fi

# Ensure all external binaries have execute permissions
echo_info "Ensuring all external binaries have execute permissions."
chmod +x "${BIN_DIR}"/*

# Optionally, add external/bin to PATH
echo_info "Dependencies have been successfully installed and configured."
echo_info "You can add the 'external/bin' directory to your PATH for easier access to these tools."
echo_info "Run the following command in your terminal, or add it to your ~/.bashrc or ~/.bash_profile:"
echo_info "export PATH=\"${BIN_DIR}:\$PATH\""
echo_info "Then, reload your shell configuration with:"
echo_info "source ~/.bashrc"

echo_info "Setup completed. Please ensure all scripts have been correctly updated and rerun your pipeline."
