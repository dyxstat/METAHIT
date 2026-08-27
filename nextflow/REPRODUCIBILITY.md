# METAHICT native DSL2 reproducibility record

This is the canonical human-readable record for the native Nextflow workflow.
It is intended to keep the README, workflow code, environment locks, container
recipes, Methods text, and Supplementary Information aligned.

## Workflow implementation

`main_dsl2.nf` is the entry point and
`modules/local/metahict_modules.nf` contains the named DSL2 processes.  Each
process calls the corresponding existing METAHICT module command through
`metahict.py`; Nextflow does not replace the underlying analysis code.

| Process | METAHICT implementation | CPUs | Memory | Time |
| --- | --- | ---: | ---: | ---: |
| `PREPROCESSING` | Module 1, once for SG and once for Hi-C | 8 | 16 GB | 12 h |
| `ASSEMBLY` | Module 2 | 32 | 128 GB | 48 h |
| `ALIGNMENT` | Module 3 | 24 | 64 GB | 24 h |
| `COVERAGE` | Module 4 | 24 | 64 GB | 24 h |
| `CONTACT` | Module 5 | 16 | 48 GB | 24 h |
| `BINNING` | Module 6 | 32 | 128 GB | 72 h |
| `REASSEMBLY` | Module 7 | 32 | 128 GB | 72 h |
| `ANNOTATION` | Module 9 | 16 | 64 GB | 48 h |
| `SCAFFOLDING` | Module 8 | 32 | 128 GB | 72 h |
| `MGE_ALIGNMENT` | Module 3 on reassembled contigs | 24 | 64 GB | 24 h |
| `MGE_CONTACT` | Module 5 on the MGE alignment | 16 | 48 GB | 24 h |
| `MGE` | Module 10 | 16 | 64 GB | 72 h |

The MGE branch consumes Module 7 `combined_contigs.fa`, realigns Hi-C reads,
and reconstructs its contact matrix before Module 10.  It does not rerun
reassembly or reuse the original-assembly contact matrix.

## Defaults and parameters

The native workflow preserves the analytical defaults of the numbered module
scripts.  It does not silently substitute a second set of defaults.  The
human-readable mirror is `assets/default_params.yaml`; the runtime authority
is the invoked module script and its command-line help.

The Nextflow workflow fixes the following stage-level choices:

- Module 1 runs without deduplication for shotgun reads and with `--dedup` for
  Hi-C reads.
- Module 2 is invoked with `--megahit`.
- Module 5 and `MGE_CONTACT` are invoked with `normcc`.
- The default restriction-enzyme string is `Sau3AI,MluCI` when the sample sheet
  does not supply an `enzyme` value.
- Module 6 is invoked with `--no-spades` because its input is the Module 2
  assembly.

Per-sample changes are supplied through the documented `*_extra_args` sample
sheet columns.  Required input paths and output directories are supplied by
Nextflow, not by `default_params.yaml`.

## Software and databases

`installation/locks/linux-64/` is the authoritative exact Linux dependency
record.  It contains one explicit artifact lock for each runtime environment:
`metahict_env`, `checkm2`, `gtdbtk-2.4.0`, `genomad`, and `ccfind_env`.
`installation/pip-requirements.txt` pins the METAHICT-compatible bin3C,
ImputeCC, and MetaCC code to immutable Git commits.

`THIRD_PARTY.md` records the upstream source, licence, version/build or
commit, and citation for direct scientific dependencies.  The Dockerfile and
Apptainer definition create the same named environments from these locks.

External databases are never baked into the container image.  The distributed
database scripts select GTDB release 220.  geNomad,
CheckM, and CheckM2 database versions are recorded by their official download
tools in the installed database directories.

External database locations are ordinary Nextflow parameters:
`--checkm_db`, `--checkm2_db`, `--gtdbtk_db`, and `--genomad_db`.
The supplied defaults resolve to `<project>/databases/`; a
custom location is supplied by overriding the relevant parameter.  The
workflow forwards each value to the modules that need it, including CheckM2,
so users do not export container-specific environment variables.

## Execution and provenance

Users download the METAHICT repository and choose one execution profile:

- `conda`: uses the exact project-local environment bundle after
  `installation/run_setup_in_venv.sh`.
- `docker`: runs the declared container image.
- `apptainer`: runs the declared SIF/container image on HPC systems.

Docker and Apptainer are alternatives, not additional installations.  The
container image supplies software only; the downloaded repository supplies the
Nextflow workflow and METAHICT scripts.  The release-image recipe and future
publication procedure are documented in `containers/PUBLISHING.md`.

Every execution writes the Nextflow trace, timeline, report, directed acyclic
graph, and parameter snapshot under the configured report directory.  The
5% human-gut test dataset is documented in `../docs/test_dataset.md`.

## Manuscript and Supplementary Information text

Use [MANUSCRIPT_TEXT.md](MANUSCRIPT_TEXT.md) as the source text for the
Methods and Code Availability sections.  Before submission, replace its
explicit final-release placeholders with the tested commit, public image
digest, GitHub Release link, and permanent archive DOI.  No release identifier
is invented in this working tree.
