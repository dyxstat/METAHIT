#!/bin/bash

# run_setup_in_venv.sh
# Sets up and activates a minimal Conda environment to run setup.sh.

set -e

ENV_NAME="metahict_venv"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOCK_FILE="${SCRIPT_DIR}/locks/linux-64/${ENV_NAME}.explicit.txt"
ENV_PREFIX="${PROJECT_ROOT}/conda_envs/${ENV_NAME}"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "[ERROR] The distributed release locks target Linux only. Use a tested platform-specific lock." >&2
    exit 1
fi

if [[ ! -s "$LOCK_FILE" ]]; then
    echo "[ERROR] Missing bootstrap lock: $LOCK_FILE" >&2
    exit 1
fi

echo "[INFO] Checking if Conda environment '$ENV_NAME' exists..."
if [[ -d "$ENV_PREFIX" ]]; then
    installed_sha256="$(conda list --explicit -p "$ENV_PREFIX" | sha256sum | awk '{print $1}')"
    expected_sha256="$(sha256sum "$LOCK_FILE" | awk '{print $1}')"
    if [[ "$installed_sha256" != "$expected_sha256" ]]; then
        echo "[ERROR] Existing '$ENV_NAME' does not match the release lock. Remove it, then rerun this installer." >&2
        exit 1
    fi
    echo "[INFO] Environment '$ENV_NAME' matches the release lock. Skipping creation."
else
    echo "[INFO] Creating new Conda environment: $ENV_NAME"
    conda create -p "$ENV_PREFIX" -y --file "$LOCK_FILE"
fi

# Locate conda.sh dynamically
if [ -z "$CONDA_EXE" ]; then
    echo "[ERROR] Conda not found in PATH. Please install Conda or add it to your PATH."
    exit 1
fi

# Derive conda.sh path automatically
CONDA_SH="$(dirname $(dirname $CONDA_EXE))/etc/profile.d/conda.sh"

if [ ! -f "$CONDA_SH" ]; then
    echo "[ERROR] Could not locate conda.sh automatically at $CONDA_SH."
    echo "Please check your Conda installation."
    exit 1
fi

# Source conda.sh and activate environment
echo "[INFO] Activating Conda environment '$ENV_NAME'..."
source "$CONDA_SH"
conda activate "$ENV_PREFIX"


# Run the setup script
echo "[INFO] Running setup.sh inside '$ENV_NAME'..."
(cd "$SCRIPT_DIR" && bash setup.sh)

echo "[INFO] setup.sh completed successfully in environment '$ENV_NAME'."
