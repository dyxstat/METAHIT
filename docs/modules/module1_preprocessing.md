# Module 1: preprocessing

Preprocessing cleans paired-end shotgun and Hi-C reads with BBTools/BBDuk and produces FastQC reports before and after trimming unless those reports are skipped.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- |
| `sg1`, `sg2` | Nextflow sample sheet | Shotgun paired-end FASTQ files |
| `hic1`, `hic2` | Nextflow sample sheet | Hi-C paired-end FASTQ files |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `-1` | module wrapper | Forward FASTQ file |
| `-2` | module wrapper | Reverse FASTQ file |
| `-o`, `--output` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--dedup` | shotgun: `false`; Hi-C: `true` | Enable duplicate removal |
| `--prefix` | derived from input name | Custom prefix for output files |
| `--minlen` | `50` | Minimum read length retained after trimming |
| `--trimq` | `10` | Quality-trimming threshold |
| `--qtrim` | `r` | BBDuk quality-trimming direction |
| `--ftl` | `10` | Number of bases trimmed from the left end |
| `--xmx` | 80% of available memory | Java memory for BBDuk/BBMap tools |
| `--ftm` | `5` | BBDuk modulo trimming value |
| `--ktrim` | `r` | Adapter-trimming direction |
| `--k` | `23` | Adapter k-mer size |
| `--mink` | `11` | Minimum adapter k-mer size |
| `--hdist` | `1` | Adapter k-mer Hamming distance |
| `--adapter-ref` | - | Adapter reference FASTA; default uses the bundled adapter reference |
| `--skip-pre-qc-report` | `false` | Skip FastQC report for input reads |
| `--skip-post-qc-report` | `false` | Skip FastQC report for cleaned reads |

## Outputs

| Output | Description |
| --- | --- |
| `1_preprocessing/sg/` | Cleaned shotgun reads |
| `1_preprocessing/hic/` | Cleaned Hi-C reads |
| `final_*_1.fastq.gz`, `final_*_2.fastq.gz` | Cleaned paired-end FASTQ files |
| FastQC reports | Pre- and post-cleaning QC reports unless skipped |

## Nextflow command

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module1 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/module1" \
  --report_dir "$PWD/results/module1/nextflow_reports" \
  -work-dir "$PWD/results/module1/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py preprocessing \
  -p "$PWD" \
  -1 example_data/sg_R1.fastq.gz \
  -2 example_data/sg_R2.fastq.gz \
  -o results/manual/module1_sg \
  -t 80
```

Run the wrapper a second time for the matched Hi-C reads:

```bash
python metahict.py preprocessing \
  -p "$PWD" \
  -1 example_data/hic_R1.fastq.gz \
  -2 example_data/hic_R2.fastq.gz \
  -o results/manual/module1_hic \
  -t 80 \
  --dedup
```
