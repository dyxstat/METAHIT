# Command-by-command tutorial

This tutorial runs one sample through the complete METAHICT workflow. The
shotgun and Hi-C libraries must represent the same biological sample.

## Commands at a glance

```bash
git clone https://github.com/dyxstat/METAHICT.git
cd METAHICT

./metahict doctor
./metahict install
./metahict test workflow

./metahict database all
./metahict doctor --runtime --databases

./metahict config
./metahict samplesheet --help
```

After creating `samplesheet.csv`, run the analysis:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results
```

The following sections explain each step.

## 1. Check the server and input data

METAHICT 1.2.0 supports and tests on 64-bit x86 Linux systems. The host must
provide Conda, `curl`, `tar`, and Git. Python, Nextflow, and the scientific
programs are installed into project-local environments by METAHICT.

Prepare these sample facts:

| Field | Meaning |
| --- | --- |
| `sample` | Short unique identifier without spaces |
| `sg1`, `sg2` | Paired short-read shotgun FASTQs, or one long-read file in `sg1` with empty `sg2` |
| `hic1`, `hic2` | Paired metagenomic Hi-C FASTQ files |
| `enzyme` | Restriction enzyme or comma-separated enzyme combination used for the Hi-C library |
| `long_read_type` | Empty for short reads; metaFlye type for a long-read `sg1` file |

The enzyme must come from the library-preparation record; there is no reliable default.

## 2. Download METAHICT

```bash
git clone https://github.com/dyxstat/METAHICT.git
cd METAHICT
```

All following relative commands assume that the current directory is the
METAHICT checkout.

## 3. Install and test METAHICT

```bash
./metahict doctor
./metahict install
./metahict test workflow
```

Expected outcome:

- `doctor` confirms the supported host and distributed lock checksums;
- `install` creates and verifies the project-local environments under
  `conda_envs/`; and
- `test workflow` runs the core workflow and standalone scaffolding entry with
  dependency-free stub processes.

The workflow test does not need the reference databases and does not perform a
biological analysis. Continue only after `install` and `test workflow` pass.

## 4. Install the databases

```bash
./metahict database all
./metahict doctor --runtime --databases
```

The default installation places CheckM, CheckM2, GTDB-Tk, and geNomad data
under `databases/`. Existing or shared installations can be reused as
described in [Installation and databases](installation.md).

## 5. Create the YAML configuration

```bash
./metahict config
```

This creates `metahict_configuration.yaml`. The file has two main sections:

- `resources`: threads and memory for each workflow stage;
- `modules`: scientific methods and their parameters.

For a first run, keep the scientific defaults. Compare the resource settings
with the server before starting. See [Configuration](configuration.md) for the
complete parameter reference.

## 6. Create the samplesheet

```bash
./metahict samplesheet \
  --sample sample_01 \
  --sg-r1 /data/sample_01/shotgun_R1.fastq.gz \
  --sg-r2 /data/sample_01/shotgun_R2.fastq.gz \
  --hic-r1 /data/sample_01/hic_R1.fastq.gz \
  --hic-r2 /data/sample_01/hic_R2.fastq.gz \
  --enzyme Sau3AI,MluCI \
  --output samplesheet.csv
```

Replace every example value. The generated CSV stores absolute FASTQ paths.
For another sample, repeat the command with `--append` and the same output
file.

For a long-read shotgun sample, omit `--sg-r2` and identify the metaFlye input
type:

```bash
./metahict samplesheet \
  --sample long_sample_01 \
  --sg-r1 /data/long_sample_01/shotgun_long_reads.fastq.gz \
  --long-read-type nano-hq \
  --hic-r1 /data/long_sample_01/hic_R1.fastq.gz \
  --hic-r2 /data/long_sample_01/hic_R2.fastq.gz \
  --enzyme Sau3AI,MluCI \
  --output samplesheet.csv
```

Valid types are `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`,
`nano-corr`, and `nano-hq`.

Check the result:

```bash
column -s, -t samplesheet.csv
```

## 7. Run the complete workflow

Start the complete workflow:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results
```

A successful invocation ends with `[PASS] METAHICT workflow completed`.
Scaffolding is not run automatically. It is an optional per-MAG analysis; use
`./metahict run --entry-module scaffolding --help` after selecting a MAG.

## 8. Check the main results

For `sample_01`, confirm the main assembly and MAG products:

```bash
test -s results/sample_01/2_assembly/final_assembly.fasta
test -d results/sample_01/6_binning/metahict/final_bins
test -s results/sample_01/6_binning/metahict/final_bins_quality.tsv
```

No output means that these paths exist. The numbered directories contain the
stage results; `results/nextflow_reports/` contains logs and provenance.
Review the quality tables before biological interpretation.

For a paired short-read sample, also check
`7_reassembly/reassembled_bins/` and `7_reassembly/reassembled_bins_quality.tsv`.
Long-read samples skip stage 7 and use the final MAGs from stage 6 in the
downstream stages.

## Resume after a failure

Correct the reported problem, keep `results/nextflow_work/`, and repeat the
same command with `--resume`:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results \
  --resume
```

Start diagnosis with [Troubleshooting](troubleshooting.md).

## Run from outside the checkout

Use the absolute path to the launcher and the analysis files:

```bash
/path/to/METAHICT/metahict run \
  --samplesheet /path/to/analysis/samplesheet.csv \
  --config /path/to/analysis/metahict_configuration.yaml \
  --outdir /path/to/analysis/results
```

The launcher locates the workflow, environments, and modules relative to its
own installation. Relative analysis paths are resolved from the current
directory.

## Next steps

- Run the bundled real-data test: [Testing METAHICT](test_dataset.md).
- Change methods or resources: [Configuration](configuration.md).
- Run one stage with existing results: [Module reference](modules/README.md).
- Interpret logs and provenance: [Logging](logging.md).
