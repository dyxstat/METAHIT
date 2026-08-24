# Example dataset

The METAHICT example dataset is a 5% subset of the human-gut benchmark dataset. The repository provides metadata, checksums, and a Nextflow sample sheet for this dataset; the paired-end shotgun and metagenomic Hi-C FASTQ files are archived externally and must be downloaded from Zenodo before running the example workflow.

The example FASTQ files are archived at:

```text
https://doi.org/10.5281/zenodo.21695166
```

Place the files under `example_data/` with these names:

```text
example_data/sg_R1.fastq.gz
example_data/sg_R2.fastq.gz
example_data/hic_R1.fastq.gz
example_data/hic_R2.fastq.gz
```

The repository includes:

| File | Purpose |
| --- | --- |
| `example_data/README.md` | Example-data description |
| `example_data/manifest.tsv` | Read counts for the example FASTQ files |
| `example_data/MD5SUMS.txt` | Checksums for downloaded example files |
| `nextflow/assets/example_data_samplesheet.csv` | Nextflow sample sheet for the example run |

## Run the smoke test

```bash
bash nextflow/ci/run_smoke_ci.sh
```

## Run the example dataset

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/example_data" \
  --report_dir "$PWD/results/example_data/nextflow_reports" \
  -work-dir "$PWD/results/example_data/nextflow_work" \
  -ansi-log false
```

## Expected output layout

For a sample named `example_data`, the output is written under:

```text
results/example_data/example_data/
```

Main subdirectories:

| Directory | Content |
| --- | --- |
| `1_preprocessing/sg/` | Cleaned shotgun reads and QC reports |
| `1_preprocessing/hic/` | Cleaned Hi-C reads and QC reports |
| `2_assembly/` | Assembled contigs |
| `3_alignment/` | Hi-C alignments and alignment diagnostics |
| `4_coverage/` | Shotgun coverage table |
| `5_contact/` | Raw and normalized Hi-C contact matrices |
| `6_binning/` | Hi-C binner outputs and consolidated MAGs |
| `7_reassembly/` | Reassembled bins and combined contig set |
| `8_scaffolding/` | Hi-C-guided scaffold outputs and contact heatmaps |
| `9_annotation/` | GTDB-Tk annotations |
| `10_MGE/` | MGE calls, candidate MGE-MAG associations, and sequence topology |
