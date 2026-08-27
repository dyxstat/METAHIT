# Test dataset

The METAHICT test dataset is a 5% subset of the human-gut benchmark dataset. The repository provides the Nextflow sample sheet and expected-output checks for this dataset; the paired-end shotgun and metagenomic Hi-C FASTQ files, metadata, and checksums are archived externally and must be downloaded from Zenodo before running the test workflow.

The test FASTQ files are archived at:

```text
https://doi.org/10.5281/zenodo.21695166
```

Place the files under `test_data/` with these names:

```text
test_data/sg_R1.fastq.gz
test_data/sg_R2.fastq.gz
test_data/hic_R1.fastq.gz
test_data/hic_R2.fastq.gz
```

The repository includes:

| File | Purpose |
| --- | --- |
| `nextflow/assets/test_data_samplesheet.csv` | Nextflow sample sheet for the test run |
| `nextflow/tests/expected/test_data_outputs.tsv` | Expected-output checks for the test run |

The Zenodo archive includes `README.md`, `manifest.tsv`, `MD5SUMS.txt`, and the four FASTQ files.

## Run the smoke test

```bash
bash nextflow/ci/run_smoke_ci.sh
```

## Run the test dataset

```bash
nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --out_root "$PWD/results/test_data" \
  --report_dir "$PWD/results/test_data/nextflow_reports" \
  -work-dir "$PWD/results/test_data/nextflow_work" \
  -ansi-log false
```

## Expected output layout

For a sample named `test_data`, the output is written under:

```text
results/test_data/test_data/
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
