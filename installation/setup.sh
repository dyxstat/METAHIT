#!/bin/bash

# setup.sh
# This script sets up the necessary dependencies for the METAHICT pipeline.
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

# Resolve all paths from this script, rather than from the caller's directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOCK_DIR="${SCRIPT_DIR}/locks/linux-64"
EXTERNAL_DIR="${PROJECT_ROOT}/external"
BIN_DIR="${EXTERNAL_DIR}/bin"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo_error "The distributed release locks target Linux only. Use a tested platform-specific lock."
    exit 1
fi

create_locked_environment() {
    local env_name="$1"
    local lock_file="${LOCK_DIR}/${env_name}.explicit.txt"
    if [[ ! -s "$lock_file" ]]; then
        echo_error "Missing lock file: $lock_file"
        exit 1
    fi
    local env_prefix="${PROJECT_ROOT}/conda_envs/${env_name}"
    if [[ -d "$env_prefix" ]]; then
        local installed_sha256
        local expected_sha256
        installed_sha256="$(conda list --explicit -p "$env_prefix" | sha256sum | awk '{print $1}')"
        expected_sha256="$(sha256sum "$lock_file" | awk '{print $1}')"
        if [[ "$installed_sha256" != "$expected_sha256" ]]; then
            echo_error "Existing '${env_name}' does not match ${lock_file##*/}. Remove that environment, then rerun this installer."
            exit 1
        fi
        echo_info "Environment '${env_name}' matches ${lock_file##*/}. Skipping creation."
    else
        echo_info "Creating '${env_name}' from ${lock_file##*/}"
        conda create -y -p "$env_prefix" --file "$lock_file"
    fi
}

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

function install_bbtools() {
    BBTOOLS_VERSION="39.10"
    BBTOOLS_TARBALL="Bbmap_${BBTOOLS_VERSION}.tar.gz"
    BBTOOLS_URL="https://sourceforge.net/projects/bbmap/files/BBMap_${BBTOOLS_VERSION}.tar.gz/download"
    BBTOOLS_SHA256="ab5dfc0bbaa5be338596aec3558c7a7c891e8d8b186e9bd671552466215b9b15"
    if [ ! -f "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" ]; then
        echo_info "Downloading BBTools ${BBTOOLS_VERSION}..."
        wget -O "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" "${BBTOOLS_URL}"
    else
        echo_info "BBTools tarball already exists, skipping download."
    fi

    printf '%s  %s\n' "$BBTOOLS_SHA256" "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" | sha256sum -c -

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



create_locked_environment "gtdbtk-2.4.0"
create_locked_environment "metahict_env"
create_locked_environment "checkm2"
create_locked_environment "genomad"
create_locked_environment "checkv_env"

echo_info "Installing pinned Pip dependencies for 'metahict_env'..."
conda run -p "${PROJECT_ROOT}/conda_envs/metahict_env" python -m pip install --no-deps --upgrade --force-reinstall \
    -r "${SCRIPT_DIR}/pip-requirements.txt"

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
