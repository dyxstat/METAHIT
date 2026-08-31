# METAHICT continuous integration

This directory contains the source, interface, and workflow-stub checks used
to validate METAHICT.

## Local CI Smoke Test

Run the local smoke test from the METAHICT repository root:

```bash
./metahict test workflow
```

The smoke test uses Nextflow stub mode. It checks workflow parsing, module
ordering, profile loading, report generation, and sample-sheet handling without
requiring sequencing data or large reference databases. Tiny valid paired
short-read and single-file long-read FASTQs and their samplesheet are created
in the temporary test directory.

## GitHub Actions

The active workflow is stored at:

```text
.github/workflows/metahict-smoke.yml
```

It runs the source-policy, unit, syntax, lock-checksum, and Nextflow stub tests.
Full biological validation uses the bundled example dataset on a Linux system
with the locked environments and reference databases installed.
