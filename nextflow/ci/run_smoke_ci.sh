#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEXTFLOW_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_DIR="$(cd "${NEXTFLOW_DIR}/.." && pwd)"

cd "$PROJECT_DIR"

export PATH="${NEXTFLOW_DIR}/bin:${PROJECT_DIR}/conda_envs/metahict_venv/bin:${PATH}"
export NXF_HOME="${NEXTFLOW_DIR}/.nextflow"

nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --out_root "${NEXTFLOW_DIR}/test_runs/ci_smoke/results" \
  --report_dir "${NEXTFLOW_DIR}/reports/ci_smoke" \
  -work-dir "${NEXTFLOW_DIR}/work/ci_smoke" \
  --threads 2 \
  --clean true \
  --chain true \
  -stub-run \
  -ansi-log false

python3 nextflow/bin/check_expected_outputs.py \
  --root "${NEXTFLOW_DIR}/test_runs/ci_smoke/results" \
  --manifest nextflow/tests/expected/native_dsl2_stub_outputs.tsv
