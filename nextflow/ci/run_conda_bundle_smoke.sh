#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEXTFLOW_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_DIR="$(cd "${NEXTFLOW_DIR}/.." && pwd)"

cd "$PROJECT_DIR"

export PATH="${NEXTFLOW_DIR}/bin:${PROJECT_DIR}/conda_envs/metahict_venv/bin:${PATH}"
export NXF_HOME="${NEXTFLOW_DIR}/.nextflow"

bash installation/verify_env_bundle.sh --project "$PROJECT_DIR"

nextflow run nextflow/main_dsl2.nf \
  -profile conda \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --out_root "${NEXTFLOW_DIR}/test_runs/conda_bundle_smoke/results" \
  --report_dir "${NEXTFLOW_DIR}/reports/conda_bundle_smoke" \
  -work-dir "${NEXTFLOW_DIR}/work/conda_bundle_smoke" \
  -stub-run \
  -ansi-log false

python3 nextflow/bin/check_expected_outputs.py \
  --root "${NEXTFLOW_DIR}/test_runs/conda_bundle_smoke/results" \
  --manifest nextflow/tests/expected/native_dsl2_stub_outputs.tsv
