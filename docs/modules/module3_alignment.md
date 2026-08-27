# Module 3: alignment

Alignment maps cleaned Hi-C reads to assembled contigs with BWA-MEM and filters alignments for downstream contact construction.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--assembly_dir` | Nextflow selected-module run | Assembly output directory containing the contig FASTA |
| `--hic_preprocessing_dir` | Nextflow selected-module run | Hi-C preprocessing output directory |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `-r`, `--reference` | module wrapper | Assembly FASTA used as the alignment reference |
| `-1` | module wrapper | Forward preprocessed Hi-C read file |
| `-2` | module wrapper | Reverse preprocessed Hi-C read file |
| `-o`, `--output` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--bwa-options` | `-5SP` | BWA-MEM options |
| `--samtools-filter` | `-F 0x900` | SAM flag filter |
| `--mapq` | `30` | Minimum mapping quality |
| `--min-match-len` | `30` | Minimum nucleotide match length |
| `--min-intra-dist` | `10000` | Minimum intra-contig distance used in metrics |
| `--sort-memory` | `1G` | Memory per SAMtools sort thread |
| `--tmp-dir` | system temporary directory | Temporary directory root |
| `--keep-sam` | `false` | Keep intermediate SAM file |
| `--skip-metrics` | `false` | Skip alignment diagnostics |

## Outputs

| Output | Description |
| --- | --- |
| `3_alignment/sorted_map.bam` | Filtered, sorted Hi-C alignment BAM |
| alignment metrics | Alignment-derived QC metrics, including 3D ratio unless skipped |

## Nextflow command

```bash
nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module3 \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --assembly_dir "$PWD/results/test_data/test_data/2_assembly" \
  --hic_preprocessing_dir "$PWD/results/test_data/test_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module3" \
  --report_dir "$PWD/results/module3/nextflow_reports" \
  -work-dir "$PWD/results/module3/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py alignment \
  -p "$PWD" \
  -r results/manual/module2/final_assembly.fasta \
  -1 results/manual/module1_hic/final_*_1.fastq.gz \
  -2 results/manual/module1_hic/final_*_2.fastq.gz \
  -o results/manual/module3 \
  -t 80
```
