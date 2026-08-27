# Module 2: assembly

Assembly constructs contigs from cleaned shotgun reads. MEGAHIT is the default short-read assembler; metaSPAdes and metaFlye are also exposed by the wrapper.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--sg_preprocessing_dir` | Nextflow selected-module run | Shotgun preprocessing output directory |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `-1` | module wrapper | Forward preprocessed read file |
| `-2` | module wrapper | Reverse preprocessed read file |
| `-o`, `--output` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `assembler` | `megahit` | Assembler selected by the Nextflow defaults |
| `--megahit` | enabled by default | Use MEGAHIT |
| `--metaspades` | disabled by default | Use metaSPAdes |
| `--metaflye` | disabled by default | Use metaFlye |
| `--memory` | 80% of available memory | Memory limit passed to the assembler where applicable |
| `--min-len` | `1000` | Minimum contig length retained |
| `--k-min` | `21` | MEGAHIT minimum k-mer size |
| `--k-max` | `141` | MEGAHIT maximum k-mer size |
| `--k-step` | `12` | MEGAHIT k-mer step size |
| `--merge-level` | `20,0.95` | MEGAHIT merge-level setting |
| `--k-list` | `21,33,55` | metaSPAdes k-mer list |
| `--flye-method` | `--nano-raw` | metaFlye read-type argument |
| `--tmp-dir` | system temporary directory | Temporary directory root |
| `--skip-quast` | `false` | Skip QUAST report |
| `--keep-temp` | `false` | Keep temporary assembler files |

## Outputs

| Output | Description |
| --- | --- |
| `2_assembly/final_assembly.fasta` | Assembly FASTA used by downstream modules |
| assembler output files | Native assembler outputs |
| QUAST report | Assembly QC report unless skipped |

## Nextflow command

```bash
nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module2 \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --sg_preprocessing_dir "$PWD/results/test_data/test_data/1_preprocessing/sg" \
  --out_root "$PWD/results/module2" \
  --report_dir "$PWD/results/module2/nextflow_reports" \
  -work-dir "$PWD/results/module2/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py assembly \
  -p "$PWD" \
  -1 results/manual/module1_sg/final_*_1.fastq.gz \
  -2 results/manual/module1_sg/final_*_2.fastq.gz \
  -o results/manual/module2 \
  -t 80 \
  --megahit
```
