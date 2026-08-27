# Module 7: reassembly

Reassembly selects short-insert, shotgun-like intra-contig Hi-C read pairs using an EM model, combines them with shotgun reads, and performs per-bin reassembly.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--binning_dir` | Nextflow selected-module run | Binning output directory with consolidated bins |
| `--assembly_dir` | Nextflow selected-module run | Original assembly directory |
| `--alignment_dir` | Nextflow selected-module run | Hi-C alignment directory |
| `--sg_preprocessing_dir` | Nextflow selected-module run | Cleaned shotgun reads |
| `--hic_preprocessing_dir` | Nextflow selected-module run | Cleaned Hi-C reads |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `--bin` | module wrapper | Directory containing input bins |
| `--assembly` | module wrapper | Original assembly FASTA |
| `--hic1`, `--hic2` | module wrapper | Preprocessed Hi-C paired-end reads |
| `--sg1`, `--sg2` | module wrapper | Preprocessed shotgun paired-end reads |
| `--bam` | module wrapper | Hi-C alignment BAM |
| `--outdir` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--cutoff-quantile` | `0.95` | Quantile used to define the short-insert cutoff |
| `--top-k` | `100` | Number of longest contigs used for EM fitting |
| `--min-mapq` | `30` | Minimum MAPQ for pair-level insert-size extraction |
| `--min-match-len` | `30` | Minimum match length |
| `--exclude-duplicates` | `false` | Exclude duplicate-marked alignments |
| `--write-nonselected-hic` | `false` | Write Hi-C pairs not selected by the EM cutoff |
| `--min-contig-len` | `500` | Minimum contig length retained after reassembly |
| `--strict-cut-off` | `2` | Strict mismatch cutoff for read recruitment |
| `--permissive-cut-off` | `5` | Permissive mismatch cutoff for read recruitment |
| `--contamination-penalty` | `5` | Penalty used to select representative reassembly |
| `--spades-mode` | `careful` | SPAdes mode; use `none` to disable careful mode |
| `--spades-phred-offset` | - | SPAdes PHRED offset setting |
| `--spades-extra-args` | - | Extra SPAdes arguments |
| `--memory` | 80% of available memory | Memory limit |
| `--skip-residual-assembly` | `false` | Skip assembly of residual reads |
| `--skip-checkm2` | `false` | Skip CheckM2 quality evaluation |
| `--checkm2_db` | `databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` | CheckM2 database path |
| `--tmp-dir` | system temporary directory | Temporary directory root |
| `--keep-temp` | `false` | Keep successful SPAdes and CheckM2 temporary directories |

## Outputs

| Output | Description |
| --- | --- |
| `7_reassembly/reassembly/reassembled_bins/` | Reassembled MAG FASTA files |
| `7_reassembly/reassembly/reassembled_bin_name_map.tsv` | Mapping from selected assembly names to simplified bin names |
| `7_reassembly/reassembly/combined_contigs.fa` | Combined binned and residual contigs for downstream modules |
| EM-selection outputs | Short-insert Hi-C read selection diagnostics |

## Nextflow command

```bash
nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module7 \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --binning_dir "$PWD/results/test_data/test_data/6_binning" \
  --assembly_dir "$PWD/results/test_data/test_data/2_assembly" \
  --alignment_dir "$PWD/results/test_data/test_data/3_alignment" \
  --sg_preprocessing_dir "$PWD/results/test_data/test_data/1_preprocessing/sg" \
  --hic_preprocessing_dir "$PWD/results/test_data/test_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module7" \
  --report_dir "$PWD/results/module7/nextflow_reports" \
  -work-dir "$PWD/results/module7/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py reassembly \
  -p "$PWD/modules" \
  --bin results/manual/module6/metahict/metahict_50_10_bins \
  --assembly results/manual/module2/final_assembly.fasta \
  --hic1 results/manual/module1_hic/final_*_1.fastq.gz \
  --hic2 results/manual/module1_hic/final_*_2.fastq.gz \
  --sg1 results/manual/module1_sg/final_*_1.fastq.gz \
  --sg2 results/manual/module1_sg/final_*_2.fastq.gz \
  --bam results/manual/module3/sorted_map.bam \
  --outdir results/manual/module7 \
  -t 80 \
  --checkm2_db "$PWD/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
```
