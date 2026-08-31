# METAHICT example dataset

This directory contains a small paired-read dataset for checking a complete
METAHICT installation with real scientific programs. It is intended for
functional testing and tutorials, not for estimating biological performance.

The shotgun reads originate from SRR6131123 and the Hi-C reads from
SRR6131122. 

| File | Library | Mate | Read pairs |
| --- | --- | --- | ---: |
| `sg_R1.fastq.gz` | Shotgun | R1 | 199,981 |
| `sg_R2.fastq.gz` | Shotgun | R2 | 199,981 |
| `hic_R1.fastq.gz` | Hi-C | R1 | 92,384 |
| `hic_R2.fastq.gz` | Hi-C | R2 | 92,384 |

The accompanying samplesheet records the restriction enzymes as
`Sau3AI,MluCI`. Run the dataset from the repository root:

```bash
./metahict run \
  --samplesheet nextflow/assets/example_dataset_samplesheet.csv \
  --config nextflow/assets/example_dataset_configuration.yaml \
  --outdir results \
  --check-outputs
```

Results are written under `results/example_dataset/`. The expected-output
check verifies the documented stage artifacts and is not a scientific
accuracy benchmark.
