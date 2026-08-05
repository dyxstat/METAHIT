# METAHICT Nextflow test assets

`expected/native_dsl2_stub_outputs.tsv` is used by the CI smoke test. It
asserts that every native DSL2 process publishes its declared result artifact,
including MGE-specific re-alignment and contact construction.

`expected/example_data_outputs.tsv` is for a real, containerized example-data
run on a system with all required reference databases. Run it with:

```bash
python nextflow/bin/check_expected_outputs.py \
  --root /path/to/metahict_results \
  --manifest nextflow/tests/expected/example_data_outputs.tsv
```

The real manifest checks non-empty biological outputs; it is not a scientific
benchmark. The bundled tiny example reads are suitable for structural stub
testing, while an end-to-end binning-to-MGE run requires sufficient sequencing
depth and configured external databases.
