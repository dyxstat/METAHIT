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

create_environment_from_yaml_if_needed() {
    local env_name="$1"
    local yaml_file="${SCRIPT_DIR}/${env_name}.yaml"
    local env_prefix="${PROJECT_ROOT}/conda_envs/${env_name}"
    if [[ -d "$env_prefix" ]]; then
        echo_info "Environment '${env_name}' already exists. Skipping creation."
        return
    fi
    if [[ ! -s "$yaml_file" ]]; then
        echo_error "Missing environment specification: $yaml_file"
        exit 1
    fi
    echo_info "Creating '${env_name}' from ${yaml_file##*/}"
    conda env create -p "$env_prefix" -f "$yaml_file"
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

function install_ccfind() {
    CCFIND_VERSION="1.4.7"
    CCFIND_COMMIT="674366b49dd31cb909c2e52834e4ec8ede8919e7"
    CCFIND_TARBALL="ccfind-${CCFIND_VERSION}.tar.gz"
    CCFIND_URL="https://github.com/yosuken/ccfind/archive/${CCFIND_COMMIT}.tar.gz"
    CCFIND_DIR="${EXTERNAL_DIR}/ccfind-${CCFIND_COMMIT}"
    CCFIND_LINK="${EXTERNAL_DIR}/ccfind"

    if [ ! -f "${EXTERNAL_DIR}/${CCFIND_TARBALL}" ]; then
        echo_info "Downloading ccfind ${CCFIND_VERSION} (${CCFIND_COMMIT})..."
        wget -O "${EXTERNAL_DIR}/${CCFIND_TARBALL}" "${CCFIND_URL}"
    else
        echo_info "ccfind tarball already exists, skipping download."
    fi

    if [ ! -d "${CCFIND_DIR}" ]; then
        echo_info "Extracting ccfind..."
        tar -xzf "${EXTERNAL_DIR}/${CCFIND_TARBALL}" -C "$EXTERNAL_DIR"
        mv "${EXTERNAL_DIR}/ccfind-${CCFIND_COMMIT}" "${CCFIND_DIR}" 2>/dev/null || true
    else
        echo_info "ccfind source directory already exists, skipping extraction."
    fi

    ln -sfn "${CCFIND_DIR}" "${CCFIND_LINK}"
    chmod +x "${CCFIND_DIR}/ccfind"
    cat > "${BIN_DIR}/ccfind" <<EOF
#!/usr/bin/env bash
exec "${CCFIND_DIR}/ccfind" "\$@"
EOF
    chmod +x "${BIN_DIR}/ccfind"
    echo_info "ccfind installed at ${CCFIND_DIR}; executable linked to ${BIN_DIR}/ccfind."
}

# Install all dependencies
install_bbtools 
install_ccfind
# Verify installations
echo_info "Verifying installations..."



create_locked_environment "gtdbtk-2.4.0"
create_locked_environment "metahict_env"
create_locked_environment "checkm2"
create_locked_environment "genomad"
if [[ -s "${LOCK_DIR}/ccfind_env.explicit.txt" ]]; then
    create_locked_environment "ccfind_env"
else
    create_environment_from_yaml_if_needed "ccfind_env"
fi

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
