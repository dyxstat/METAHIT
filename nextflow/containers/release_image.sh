#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Build, inspect, or publish the single METAHICT release image.

Usage:
  nextflow/containers/release_image.sh build <tag>
  nextflow/containers/release_image.sh inspect <tag>
  nextflow/containers/release_image.sh push <tag>

Defaults:
  repository: ghcr.io/1001shiyuan/metahict

Environment overrides:
  METAHICT_CONTAINER_REPOSITORY  destination image repository
  METAHICT_CONTAINER_TAG         image tag

`push` is the only command that uploads an image.  It requires a prior
`docker login ghcr.io` and should be used only after the local test suite has
passed.
EOF
}

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 1
fi

command="$1"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
containers_dir="${repo_root}/nextflow/containers"
repository="${METAHICT_CONTAINER_REPOSITORY:-ghcr.io/1001shiyuan/metahict}"
tag="$2"
image="${repository}:${tag}"
source_repository="https://github.com/dyxstat/METAHICT"

case "${command}" in
    build|inspect|push) ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unsupported command: ${command}" >&2; usage >&2; exit 1 ;;
esac

command -v docker >/dev/null 2>&1 || {
    echo "[ERROR] docker is required but was not found in PATH." >&2
    exit 1
}

if [[ "${command}" == "build" ]]; then
    vcs_ref="$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf 'unknown')"
    echo "[INFO] Building ${image}"
    docker build \
        --file "${containers_dir}/Dockerfile" \
        --build-arg 'METAHICT_ENVS=metahict_env checkm2 gtdbtk-2.4.0 genomad checkv_env' \
        --build-arg "METAHICT_VERSION=${tag}" \
        --build-arg "VCS_REF=${vcs_ref}" \
        --build-arg "SOURCE_REPOSITORY=${source_repository}" \
        --tag "${image}" \
        "${repo_root}"
elif [[ "${command}" == "inspect" ]]; then
    docker image inspect "${image}" >/dev/null
    echo "[INFO] Local image exists: ${image}"
    echo "[INFO] OCI metadata:"
    docker image inspect --format '{{ range $key, $value := .Config.Labels }}{{ printf "%s=%s\\n" $key $value }}{{ end }}' "${image}" | sort
else
    echo "[INFO] Publishing ${image}"
    docker push "${image}"
    echo "[INFO] Immutable digest (copy this into the release configuration):"
    docker buildx imagetools inspect "${image}" --format '{{ .Name }}@{{ .Manifest.Digest }}'
fi
