#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  nextflow/containers/build_containers.sh docker [container_key ...]
  nextflow/containers/build_containers.sh apptainer [container_key ...]
  nextflow/containers/build_containers.sh singularity [container_key ...]

container_key values are listed in nextflow/containers/module-envs.tsv.
If no container_key is provided, all unique images are built.
EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

runtime="$1"
shift

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
containers_dir="${repo_root}/nextflow/containers"
env_table="${containers_dir}/module-envs.tsv"
image_dir="${containers_dir}/images"
tag="${METAHICT_CONTAINER_TAG:-local}"
local_apptainer_bin="${containers_dir}/conda_envs/metahict_apptainer/bin"

export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${containers_dir}/apptainer_cache}"
export PATH="${local_apptainer_bin}:${PATH}"

case "$runtime" in
    docker|apptainer|singularity) ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unsupported runtime: $runtime" >&2; usage; exit 1 ;;
esac

if [[ "$runtime" == "docker" ]]; then
    command -v docker >/dev/null 2>&1 || { echo "[ERROR] docker not found in PATH." >&2; exit 1; }
else
    command -v "$runtime" >/dev/null 2>&1 || { echo "[ERROR] $runtime not found in PATH." >&2; exit 1; }
    mkdir -p "$image_dir"
fi

declare -A requested=()
if [[ $# -gt 0 ]]; then
    for key in "$@"; do
        requested["$key"]=1
    done
fi

declare -A seen_images=()
tail -n +2 "$env_table" | while IFS=$'\t' read -r key image envs; do
    if [[ ${#requested[@]} -gt 0 && -z "${requested[$key]:-}" ]]; then
        continue
    fi
    if [[ -n "${seen_images[$image]:-}" ]]; then
        continue
    fi
    seen_images["$image"]=1
    echo "[INFO] Building ${image}:${tag} for ${key} with envs: ${envs}"
    if [[ "$runtime" == "docker" ]]; then
        docker build \
            -f "${containers_dir}/Dockerfile" \
            --build-arg "METAHICT_ENVS=${envs}" \
            -t "${image}:${tag}" \
            "$repo_root"
    else
        build_def="${containers_dir}/.apptainer.${image}.${tag}.def"
        python3 "${containers_dir}/make_apptainer_def.py" \
            --template "${containers_dir}/apptainer.def" \
            --lock-dir "${repo_root}/installation/locks/linux-64" \
            --pip-requirements "${repo_root}/installation/pip-requirements.txt" \
            --out "$build_def"
        "$runtime" build \
            --force \
            --bind "${repo_root}/installation:/tmp/metahict-installation:ro" \
            --build-arg "METAHICT_ENVS=${envs}" \
            --build-arg "METAHICT_LOCK_DIR=/tmp/metahict-installation/locks/linux-64" \
            --build-arg "METAHICT_PIP_REQUIREMENTS=/tmp/metahict-installation/pip-requirements.txt" \
            "${image_dir}/${image}_${tag}.sif" \
            "$build_def"
    fi
done
