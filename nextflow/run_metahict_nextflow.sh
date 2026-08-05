#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

export NXF_HOME="${SCRIPT_DIR}/.nextflow"
export NXF_OFFLINE="${NXF_OFFLINE:-true}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SCRIPT_DIR}/containers/apptainer_cache}"
export PATH="${SCRIPT_DIR}/containers/conda_envs/metahict_apptainer/bin:${PROJECT_DIR}/conda_envs/metahict_venv/bin:${PROJECT_DIR}/conda_envs/metahict_env/bin:${PATH}"

mkdir -p "${SCRIPT_DIR}/logs"

nss_wrapper_lib="${SCRIPT_DIR}/containers/conda_envs/metahict_apptainer/lib/libnss_wrapper.so"
if ! getent passwd "$(id -u)" >/dev/null 2>&1 && [[ -s "${nss_wrapper_lib}" ]]; then
    nss_dir="${SCRIPT_DIR}/containers/nss_wrapper"
    mkdir -p "${nss_dir}"
    user_name="${USER:-metahict}"
    user_home="${HOME:-/tmp}"
    printf '%s:x:%s:%s:%s:%s:/bin/bash\n' "${user_name}" "$(id -u)" "$(id -g)" "${user_name}" "${user_home}" > "${nss_dir}/passwd"
    printf '%s:x:%s:\n' "${user_name}" "$(id -g)" > "${nss_dir}/group"
    export NSS_WRAPPER_PASSWD="${nss_dir}/passwd"
    export NSS_WRAPPER_GROUP="${nss_dir}/group"
    export LD_PRELOAD="${nss_wrapper_lib}${LD_PRELOAD:+:${LD_PRELOAD}}"
fi

exec "${SCRIPT_DIR}/bin/nextflow" -log "${SCRIPT_DIR}/logs/.nextflow.log" "$@"
