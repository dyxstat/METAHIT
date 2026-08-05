#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 --project /path/to/METAHICT" >&2
}

if [[ $# -ne 2 || "$1" != "--project" ]]; then
    usage
    exit 2
fi

project_path="$(cd "$2" && pwd)"
lock_dir="${project_path}/installation/locks/linux-64"
env_dir="${project_path}/conda_envs"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "[ERROR] METAHICT's distributed Conda locks target Linux only." >&2
    exit 1
fi

command -v conda >/dev/null 2>&1 || {
    echo "[ERROR] Conda is required for the METAHICT Conda profile." >&2
    exit 1
}

for env_name in metahict_env checkm2 gtdbtk-2.4.0 genomad checkv_env; do
    lock_file="${lock_dir}/${env_name}.explicit.txt"
    env_prefix="${env_dir}/${env_name}"
    if [[ ! -s "$lock_file" ]]; then
        echo "[ERROR] Missing release lock: $lock_file" >&2
        exit 1
    fi
    if [[ ! -d "$env_prefix" ]]; then
        echo "[ERROR] Missing required METAHICT Conda environment: $env_prefix" >&2
        echo "        Run: bash ${project_path}/installation/run_setup_in_venv.sh" >&2
        exit 1
    fi
    installed_sha256="$(conda list --explicit -p "$env_prefix" | sha256sum | awk '{print $1}')"
    expected_sha256="$(sha256sum "$lock_file" | awk '{print $1}')"
    if [[ "$installed_sha256" != "$expected_sha256" ]]; then
        echo "[ERROR] ${env_name} does not match its release lock." >&2
        echo "        Remove ${env_prefix} and rerun the METAHICT installer." >&2
        exit 1
    fi
    echo "[INFO] ${env_name} matches ${lock_file##*/}"
done
