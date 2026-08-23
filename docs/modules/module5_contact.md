# Module 5: contact generation and normalization

The contact module builds raw Hi-C contact matrices from Hi-C alignments and normalizes contacts. NormCC is the default normalization method.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- |
| `--assembly_dir` | Nextflow selected-module run | Assembly output directory containing the contig FASTA |
| `--alignment_dir` | Nextflow selected-module run | Hi-C alignment output directory |
| `enzyme` | Nextflow sample sheet | Restriction enzyme list |
| `raw`, `normcc`, `hiczin`, `bin3c`, `metator` | module wrapper | Contact command subcommand; `normcc` is the workflow default |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `--bam` | module wrapper | Hi-C alignment BAM |
| `--fasta` | module wrapper | Assembly FASTA |
| `--out` | module wrapper | Output directory |
| `--enzyme` | module wrapper | Restriction enzyme list |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `method` | `normcc` | Contact-normalization method selected by the Nextflow defaults |
| `--metacc-min-signal` | `1` | Minimum contact signal during filtering |
| `--metacc-min-len` | `1000` | Minimum contig length |
| `--metacc-min-mapq` | `30` | Minimum mapping quality |
| `--metacc-min-match` | `30` | Minimum match length |
| `--spurious-contact-percent` | `5` | Lowest percentage of normalized contacts removed as potential spurious contacts |
| `--coverage-file` | - | Optional coverage table |
| `--epsilon` | `1` | Epsilon for supported normalization methods |
| `--max-iter` | `1000` | Maximum balancing iterations |
| `--tol` | `1e-6` | Balancing convergence tolerance |

## Outputs

| Output | Description |
| --- | --- |
| `5_contact/contact/` | Contact-module output directory |
| raw contact matrix | Raw contig-by-contig contact counts |
| normalized contact matrix | Bias-corrected and spurious-contact-filtered matrix |

## Nextflow command

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module5 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --alignment_dir "$PWD/results/example_data/example_data/3_alignment" \
  --out_root "$PWD/results/module5" \
  --report_dir "$PWD/results/module5/nextflow_reports" \
  -work-dir "$PWD/results/module5/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py contact normcc \
  -p "$PWD" \
  --bam results/manual/module3/sorted_map.bam \
  --fasta results/manual/module2/final_assembly.fasta \
  --out results/manual/module5 \
  --enzyme "Sau3AI,MluCI"
```
