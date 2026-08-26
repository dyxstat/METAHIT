# Text for Methods, Code Availability, and Supplementary Information

This file is a source-text template, not a claim that a final release already
exists.  Replace the bracketed final-release fields only after the full tests,
Git tag, public container publication, and archive have been completed.

## Methods: workflow implementation

METAHICT was implemented as a native Nextflow DSL2 workflow.  The workflow
defines named processes for preprocessing of shotgun and Hi-C reads, assembly,
alignment, coverage estimation, contact construction and normalization,
binning, reassembly, annotation, scaffolding, MGE-specific realignment and
contact construction, and MGE analysis.  Each process declares explicit input
and output files, resource requirements, software environment, and published
result directory.  The established numbered METAHICT module scripts remain the
scientific implementation invoked by these processes; Nextflow provides task
orchestration, provenance, restart support, and execution profiles.

Preprocessing is run separately for shotgun and Hi-C reads, with duplicate
removal enabled only for Hi-C reads.  The default assembly is MEGAHIT and
contact normalization is NormCC.  MGE analysis uses the reassembled contigs:
Hi-C reads are realigned to the Module 7 combined contigs and a new contact
matrix is constructed before Module 10; reassembly is not repeated for the
MGE branch.  The full resource and default record is in
`nextflow/REPRODUCIBILITY.md` and `nextflow/assets/default_params.yaml`.

## Methods: software and execution

The Linux release installs exact Conda package artifacts from explicit lock
files for `metahict_env`, `checkm2`, `gtdbtk-2.4.0`, `genomad`, and `ccfind_env`.
The METAHICT-compatible bin3C, ImputeCC, and MetaCC code is installed from
separately maintained GitHub repositories at immutable commits pinned in
`installation/pip-requirements.txt`.  Sources, licences, versions, and
citations for direct dependencies are listed in `THIRD_PARTY.md`.  Docker and
Apptainer container recipes create the same locked environments.  GTDB-Tk,
geNomad, CheckM, and CheckM2 reference databases are user-configured
external data and are not embedded in the containers.  Their paths are passed
as regular Nextflow parameters and are independent of the selected execution
profile.

## Code availability

The final tested METAHICT source release is available at
https://github.com/dyxstat/METAHICT, at Git commit `[FINAL COMMIT]`, released as
Git tag `[FINAL TAG]`, and archived at `[PERMANENT DOI]`.  The tested public
container image is `[PUBLIC IMAGE REFERENCE@sha256:DIGEST]`.  The repository
contains exact Conda locks, container recipes, third-party provenance,
example data, and workflow tests.

## Supplementary Information: reproducibility record

Report the following values together: `[FINAL COMMIT]`, `[FINAL TAG]`,
`[PERMANENT DOI]`, `[PUBLIC IMAGE REFERENCE@sha256:DIGEST]`, the Nextflow
version used for testing, the selected execution profile, external database
releases/paths, and the complete example-data test result.  Do not replace
these fields with untested values or a mutable image tag.
