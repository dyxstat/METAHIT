# Installation and databases

METAHICT installs exact project-local Conda environments. The base Conda
environment is used only to launch the management command; scientific tasks
run with the locked environments selected by Nextflow.

## Supported host

METAHICT 1.2.0 supports and tests on 64-bit x86 Linux systems and requires:

- Conda;
- `curl`, `tar`, and Git;
- write access to the METAHICT checkout; and
- sufficient storage for environments, databases, results, and work files.

Check the host:

```bash
conda --version
curl --version
tar --version
git --version
```

## Download

```bash
git clone https://github.com/dyxstat/METAHICT.git
cd METAHICT
chmod +x metahict nextflow/bin/nextflow
```

## Install the runtime

```bash
./metahict doctor
./metahict install
```

`doctor` checks the platform and distributed checksums. `install` creates the
environments under `conda_envs/`. 

The installed environments are:

| Environment | Main purpose |
| --- | --- |
| `metahict_nextflow_env` | Java runtime for Nextflow |
| `metahict_env` | Core Python libraries and most bioinformatics programs |
| `metabat2` | MetaBAT2 depth summarization |
| `checkm2` | CheckM2 and DIAMOND |
| `gtdbtk-2.4.0` | GTDB-Tk classification |
| `genomad` | geNomad MGE analysis |
| `ccfind_env` | ccfind and `ssearch36` |

## Test immediately after installation

```bash
./metahict test workflow
```

The workflow test creates temporary synthetic inputs and runs the default core
workflow in stub mode. It then tests the standalone scaffolding entry with the
generated stub MAGs and checks the combined published output contract. It does
not require the reference databases. A passing result confirms the installed
workflow engine and stage connections.

Developers changing METAHICT source can additionally run
`./metahict test source`; ordinary installations do not require that test.

## Install the reference databases

The complete workflow uses four external databases:

| Database | Used for | Default location |
| --- | --- | --- |
| CheckM | Marker profiling during binning | `databases/checkm_db/` |
| CheckM2 | MAG quality evaluation | `databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` |
| GTDB-Tk release 220 | Taxonomic classification | `databases/gtdbtk_db/release220/` |
| geNomad | MGE detection and taxonomy | `databases/genomad_db/genomad_db/` |

Install all four in the default layout:

```bash
./metahict database all
./metahict doctor --runtime --databases
```

Install one database with `./metahict database checkm`, `checkm2`, `gtdbtk`,
or `genomad`.

## Shared database storage

Install under a shared root:

```bash
./metahict database all --root /shared/databases/metahict
```

Pass non-default locations to the workflow:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results \
  --checkm-db /shared/databases/metahict/checkm_db \
  --checkm2-db /shared/databases/metahict/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
  --gtdbtk-db /shared/databases/metahict/gtdbtk_db/release220 \
  --genomad-db /shared/databases/metahict/genomad_db/genomad_db
```

Selected modules require only their relevant databases. Their help pages list
the exact requirements:

```bash
./metahict run --entry-module annotation --help
```

## Reuse existing environments and databases

An existing METAHICT installation can be reused when it was created from the
same lock files. From a new checkout, link the complete environment and
database directories while keeping the original installation in place:

```bash
ln -s /path/to/existing/METAHICT/conda_envs conda_envs
ln -s /path/to/existing/METAHICT/databases databases
./metahict install
./metahict verify
./metahict doctor --runtime --databases
./metahict test workflow
```

Conda environments contain installation prefixes, so linking the complete
directory is safer than copying individual environment folders. The `install`
command verifies those environments and installs the pinned ccfind launcher in
the new checkout. Do not link the old `external/` directory because that
launcher can contain paths from the original checkout. Reuse is accepted only
when `verify` passes. If the locks differ, `install` reports the environment
that must be rebuilt.

For databases stored elsewhere, either link `databases/` or pass the four
database options shown above.

## Common installation failures

- `conda` not found: initialize Conda, open a new shell, and confirm
  `conda --version`.
- Environment lock mismatch: preserve anything needed, remove or move only the
  named project environment, and rerun `./metahict install`.
- Database missing: run `./metahict doctor --databases`, then install the
  database or pass its existing path.
- Workflow test cannot start Nextflow: confirm that
  `conda_envs/metahict_nextflow_env/bin/java` exists and rerun `verify`.

Continue with the [command-by-command tutorial](quickstart.md) or
[Testing METAHICT](test_dataset.md).
