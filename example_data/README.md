# METAHICT example data

This folder contains a 5% human-gut example dataset for running METAHICT as a
small reproducibility and installation test.

## Files

| File | Description |
| --- | --- |
| `sg_R1.fastq.gz` | Shotgun read 1, 6,270,000 read pairs from the human-gut dataset |
| `sg_R2.fastq.gz` | Shotgun read 2, 6,270,000 read pairs from the human-gut dataset |
| `hic_R1.fastq.gz` | Hi-C read 1, 4,295,000 read pairs from the human-gut dataset |
| `hic_R2.fastq.gz` | Hi-C read 2, 4,295,000 read pairs from the human-gut dataset |
| `manifest.tsv` | Read counts for each FASTQ file |
| `MD5SUMS.txt` | Checksums for the FASTQ files |

The example reads were subset as approximately 5% of the human-gut benchmark FASTQ files:
`HGSG1.fastq.gz`, `HGSG2.fastq.gz`, `HIC1.fastq.gz`, and `HIC2.fastq.gz`.

## Intended use

These files provide the raw shotgun and Hi-C paired-end inputs needed to start
the METAHICT workflow. They are intentionally smaller than the full benchmark
dataset so users can verify that the pipeline and software environments are
configured correctly.

Because this is a subset, outputs from assembly, binning, reassembly,
annotation, and MGE detection are intended for smoke testing only and should not
be interpreted as benchmark-scale biological results. Modules that require
external databases still require the corresponding METAHICT database setup.
