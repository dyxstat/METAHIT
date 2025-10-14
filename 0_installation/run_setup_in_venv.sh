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

# Determine conda.sh path 
CONDA_SH="/apps/anaconda3/2024.10-1/etc/profile.d/conda.sh"
if [ ! -f "$CONDA_SH" ]; then
    echo "[ERROR] Could not locate conda.sh at $CONDA_SH. Please adjust the script manually."
    exit 1
fi

# Source conda and activate the environment
echo "[INFO] Activating Conda environment '$ENV_NAME'..."
source "$CONDA_SH"
conda activate $ENV_NAME

# Run the setup script
echo "[INFO] Running setup.sh inside '$ENV_NAME'..."
bash setup.sh

echo "[INFO] setup.sh completed successfully in environment '$ENV_NAME'."
