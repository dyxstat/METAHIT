# Installation and databases

This page describes the local environments, database downloads, and execution profiles used by METAHICT v1.1.0.

## Local Conda installation

From the repository root:

```bash
bash installation/run_setup_in_venv.sh
```

The script creates project-local environments under `conda_envs/`:

| Environment | Main use |
| --- | --- |
| `metahict_env` | Core METAHICT wrappers and most module tools |
| `checkm2` | CheckM2 quality evaluation |
| `gtdbtk-2.4.0` | GTDB-Tk annotation |
| `genomad` | geNomad MGE discovery |
| `ccfind_env` | ccfind circularity detection |

## Reference databases

Large reference databases are not stored in the Git repository. Download them separately:

```bash
bash installation/db/checkm_db.sh databases/checkm_db
bash installation/db/checkm2_db.sh databases/checkm2_db
bash installation/db/gtdbtk_db.sh databases/gtdbtk_db/release220
bash installation/db/genomad_db.sh databases/genomad_db
```

These databases are used as follows:

| Database | Used by |
| --- | --- |
| CheckM | Hi-C binning through ImputeCC/CheckM marker profiling |
| CheckM2 | Binning, reassembly, and scaffolding quality evaluation |
| GTDB-Tk release 220 | Taxonomic annotation |
| geNomad | Viral and plasmid discovery and MGE taxonomy |

If the databases are stored outside the default `databases/` directory, pass their paths to Nextflow:

```bash
--checkm_db /path/to/checkm_db
--checkm2_db /path/to/uniref100.KO.1.dmnd
--gtdbtk_db /path/to/gtdbtk/release220
--genomad_db /path/to/genomad_db
```

## Nextflow profiles

The profiles defined in `nextflow/nextflow.config` are:

| Profile | Use case |
| --- | --- |
| `local` | Run with the project-local environments installed by `installation/run_setup_in_venv.sh` |
| `conda` | Verify and use the locked project-local Conda environment bundle |
| `docker` | Run processes with Docker |
| `apptainer` | Run processes with Apptainer, commonly used on HPC systems |
| `singularity` | Run processes with Singularity |
| `slurm` | Submit Nextflow tasks to a SLURM executor; requires site-specific cluster configuration |

For Docker, Apptainer, or Singularity runs, use `--container_image_override` to point Nextflow to the METAHICT all-tools container image.

Example Apptainer command:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile apptainer \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/example_data" \
  --report_dir "$PWD/results/example_data/nextflow_reports" \
  -work-dir "$PWD/results/example_data/nextflow_work" \
  --container_image_override "$PWD/nextflow/containers/images/metahict-all_local.sif" \
  -ansi-log false
```
