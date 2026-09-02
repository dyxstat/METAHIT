# METAHICT Nextflow test assets

`expected/workflow_stub_outputs.tsv` is used by the CI smoke test. It
asserts that every native DSL2 process publishes its declared result artifact
across the core workflow and the separate scaffolding invocation.
Its generated samplesheet contains one paired short-read sample and one
single-file long-read sample, so the test also verifies long-read routing and
the absence of shotgun preprocessing and reassembly for that row.

`expected/example_dataset_outputs.tsv` is for a real bundled-example run on a
system with all required environments and reference databases. Run it with:

```bash
python nextflow/bin/check_expected_outputs.py \
  --root results \
  --manifest nextflow/tests/expected/example_dataset_outputs.tsv
```

The real manifest checks non-empty biological outputs; it is not by itself a
scientific benchmark. After the run, compare it with the accepted baseline:

```bash
./metahict compare \
  --baseline /path/to/accepted-results/example_dataset \
  --candidate results/example_dataset \
  --report validation/scientific-regression.tsv
```

The comparison normalizes FASTA record order, applies declared numerical
tolerances to tables and sparse matrices, and compares final bin directories.
It exits nonzero on any unexplained difference.
