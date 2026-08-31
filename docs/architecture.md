# Workflow architecture

Nextflow manages the workflow graph, Python implements the scientific stages,
and a minimal POSIX-shell launcher starts the management interface.

This developer reference describes how the command-line interface, Nextflow,
Python stages, and external tools divide responsibilities.

```text
./metahict
    └── Nextflow DSL2 graph
          ├── exact environment and task resources
          ├── Python scientific stage
          │     └── checked third-party tool command
          └── published result, logs, and provenance
```

## Control plane

The control plane is implemented in Nextflow and Python. A minimal POSIX-shell
launcher obtains its Python runtime from Conda.

| Responsibility | Owner |
| --- | --- |
| End-to-end DAG and channels | `nextflow/main_dsl2.nf` |
| Process inputs, outputs, resources, environments, and stubs | `nextflow/modules/local/metahict_modules.nf` |
| Installation and exact-lock verification | `metahict_manager.py` |
| Database downloads and path validation | `metahict_manager.py` |
| Samplesheet creation | `metahict_manager.py` |
| Workflow launch, immutable logging/provenance, resume, and output checks | `metahict_manager.py`; `docs/logging.md` |
| Numbered scientific stages | `modules/*/*.py` and documented Python helpers |

Nextflow processes launch the Python stage programs directly. CI validates
this interface with `nextflow/ci/check_architecture_policy.py`.

The small top-level `metahict` file is a POSIX-shell bootstrap. It finds
Conda's Python and executes `metahict_manager.py`; scientific work begins in
the Nextflow processes.

## Process implementations

Each numbered stage is implemented in Python and is invoked directly from its
Nextflow process. External programs are called using explicit argument vectors
with checked exit status. File transformations that were formerly expressed as
shell pipelines are implemented with Python streams or explicit connected
processes. Source policy checks reject `os.system`, `shell=True`, copied
executables, and scientific shell drivers.

Every numbered directory under `modules/` has a flat layout. Its primary entry
point and helper programs are stored directly in the same directory; there are
no module-local `scripts/`, `bin_integration/`, tool-named source
subdirectories, or shell drivers. This makes the files used by a stage visible
without navigating an undocumented hierarchy.
`nextflow/ci/check_architecture_policy.py` rejects new source subdirectories in
a numbered module. Generated `__pycache__/` directories are ignored and are
not distributed as source.

## Shell boundary

The top-level `metahict` file is the only user-facing shell launcher. Nextflow
task command blocks execute under Bash, while workflow decisions and scientific
file processing remain in Nextflow and Python. Bioconda's conventional
`build.sh` and `run_test.sh` files are confined to the staging recipe.

## Change-validation rule

Because stage implementations participate in scientific computation, every
change must pass:

1. unit tests for argument and path handling;
2. the complete Nextflow stub test;
3. an end-to-end run on the bundled example dataset; and
4. comparison of scientific outputs with the accepted baseline.

`nextflow/bin/compare_scientific_outputs.py` and
`nextflow/tests/scientific_regression_manifest.tsv` provide the standard
baseline comparison. Together, these checks validate process interfaces and
scientific changes before release.

For user commands, see [Running the workflow](nextflow.md). For the inputs and
outputs of each stage, see the [Module reference](modules/README.md).
