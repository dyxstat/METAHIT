# METAHI-T Nextflow CI Templates

This directory keeps CI templates for the native METAHI-T Nextflow workflow.
They are stored under `nextflow/` so the workflow can be reviewed without
changing the root repository layout.

## Local CI Smoke Test

Run the local smoke test from the METAHI-T repository root:

```bash
bash nextflow/ci/run_smoke_ci.sh
```

The smoke test uses Nextflow stub mode. It checks workflow parsing, module
ordering, profile loading, report generation, and sample-sheet handling without
requiring large reference databases.

## GitHub Actions Template

`github-actions-smoke.yml` is a ready-to-enable GitHub Actions workflow. To make
it active in a GitHub repository, copy it to:

```text
.github/workflows/metahict-nextflow-smoke.yml
```

Full biological tests should be run on a machine with the METAHI-T databases
and container runtime available. The GitHub Actions template is intentionally a
lightweight orchestration check.
