# METAHICT Nextflow Changelog

## Unreleased workflow refinement

- Replaced the generic wrapper execution path with native DSL2 processes for
  preprocessing, assembly, alignment, coverage, contact construction, binning,
  reassembly, scaffolding, annotation, and MGE analysis.
- Made the MGE branch explicit: it realigns Hi-C reads to reassembled contigs,
  reconstructs a contact matrix, and then runs MGE analysis.
- Added declared task resources, Conda environment, container image, inputs,
  outputs, and result publication for every process.
- Added a complete DSL2 stub-output manifest to the CI smoke test.
- Added exact Linux Conda artifact locks, a multi-environment verification
  check, and reproducible installation records.
- Added container build recipes and release-image publication instructions.
- Added third-party provenance and release/archive records.
- Preserved the former generic wrapper as `main_wrapper_legacy.nf` for
  comparison; it is not called by the native entry point.
