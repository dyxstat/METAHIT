#!/bin/bash

# run_setup_in_venv.sh
# Sets up and activates a minimal Conda environment to run setup.sh.

set -e

ENV_NAME="metahit_venv"

echo "[INFO] Checking if Conda environment '$ENV_NAME' exists..."
if conda info --envs | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "[INFO] Environment '$ENV_NAME' already exists. Skipping creation."
else
    echo "[INFO] Creating new Conda environment: $ENV_NAME"
    conda create -n $ENV_NAME -y -c bioconda -c conda-forge \
        wget \
        unzip \
        openjdk \
        perl \
        git
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
conda activate "$ENV_NAME"


# Run the setup script
echo "[INFO] Running setup.sh inside '$ENV_NAME'..."
bash setup.sh

echo "[INFO] setup.sh completed successfully in environment '$ENV_NAME'."
