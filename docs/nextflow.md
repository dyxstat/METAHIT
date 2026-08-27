# Nextflow workflow usage

The main workflow entry point is:

```text
nextflow/main_dsl2.nf
```

## Full workflow

Prepare the shell from the METAHICT repository root:

```bash
export JAVA_HOME="$PWD/conda_envs/metahict_venv"
export PATH="$PWD/nextflow/bin:$PWD/conda_envs/metahict_venv/bin:$PWD/conda_envs/metahict_env/bin:$PWD/external/bin:$PATH"
export NXF_HOME="$PWD/nextflow/.nextflow"
```

```bash
nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --out_root "$PWD/results/test_data" \
  --report_dir "$PWD/results/test_data/nextflow_reports" \
  -work-dir "$PWD/results/test_data/nextflow_work" \
  -ansi-log false
```

## Sample sheet

The sample sheet is a CSV file used by Nextflow to define sample names, raw read paths, restriction enzymes, optional pre-existing module outputs, and module-specific extra arguments.

The test sample sheet is:

```text
nextflow/assets/test_data_samplesheet.csv
```

For a raw short-read run, the essential columns are:

| Column | Description |
| --- | --- |
| `sample` | Sample name |
| `sg1`, `sg2` | Shotgun paired-end FASTQ files |
| `hic1`, `hic2` | Hi-C paired-end FASTQ files |
| `enzyme` | Restriction enzyme list, for example `Sau3AI,MluCI` |

Additional columns in the distributed sample sheets demonstrate how to provide optional inputs and module-specific extra arguments.

## Selected-module execution

Use `--entry_module` to run one workflow module at a time. Valid values are:

```text
all, module1, module2, module3, module4, module5, module6, module7, module8, module9, module10
```

Each module page gives a module-specific Nextflow command:

| Entry module | Page |
| --- | --- |
| `module1` | [Preprocessing](modules/module1_preprocessing.md) |
| `module2` | [Assembly](modules/module2_assembly.md) |
| `module3` | [Alignment](modules/module3_alignment.md) |
| `module4` | [Coverage](modules/module4_coverage.md) |
| `module5` | [Contact generation and normalization](modules/module5_contact.md) |
| `module6` | [Binning and consolidation](modules/module6_binning.md) |
| `module7` | [Reassembly](modules/module7_reassembly.md) |
| `module8` | [Scaffolding](modules/module8_scaffolding.md) |
| `module9` | [Annotation](modules/module9_annotation.md) |
| `module10` | [MGE analysis](modules/module10_mge.md) |

## Reusing outputs from previous modules

Selected-module runs can receive previous outputs through direct path parameters:

| Parameter | Expected content |
| --- | --- |
| `--sg_preprocessing_dir` | Shotgun preprocessing output directory |
| `--hic_preprocessing_dir` | Hi-C preprocessing output directory |
| `--assembly_dir` | Assembly output directory |
| `--alignment_dir` | Alignment output directory |
| `--coverage_dir` | Coverage output directory |
| `--contact_dir` | Contact output directory |
| `--binning_dir` | Binning output directory |
| `--reassembly_dir` | Reassembly output directory |
| `--mge_alignment_dir` | MGE alignment output directory |
| `--mge_contact_dir` | MGE contact output directory |

## Common workflow options

| Option | Meaning |
| --- | --- |
| `-profile` | Nextflow execution profile: `local`, `conda`, `docker`, `apptainer`, `singularity`, or `slurm` |
| `--samplesheet` | CSV file describing sample inputs |
| `--entry_module` | Workflow entry point; default is `all` |
| `--out_root` | Output root directory |
| `--report_dir` | Nextflow reports and provenance directory |
| `-work-dir` | Nextflow working directory |
| `--container_image_override` | Path or image name for the METAHICT all-tools container used by Docker, Apptainer, or Singularity profiles |
| `-ansi-log false` | Plain log output, useful for saved logs |
