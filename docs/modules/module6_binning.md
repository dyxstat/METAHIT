# Module 6: binning and consolidation

Binning runs Hi-C-aware binners, generates hybrid bin sets, evaluates candidate bins, and consolidates them into a non-redundant MAG collection.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--assembly_dir` | Nextflow selected-module run | Assembly output directory containing the contig FASTA |
| `--alignment_dir` | Nextflow selected-module run | Hi-C alignment output directory |
| `enzyme` | Nextflow sample sheet | Restriction enzyme list |
| `--project_path` | module wrapper | METAHICT module path |
| `--fasta` | module wrapper | Assembly FASTA |
| `--bam` | module wrapper | Hi-C alignment BAM |
| `--output` | module wrapper | Output directory |
| `--enzyme` | module wrapper | Restriction enzyme list |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--checkm_db` | `databases/checkm_db` | CheckM database path |
| `--metacc-min-len` | `1000` | Minimum contig length for MetaCC-style filtering |
| `--metacc-min-signal` | `2` | Minimum contact signal |
| `--metacc-min-mapq` | `30` | Minimum MAPQ for MetaCC-style filtering |
| `--metacc-min-match` | `30` | Minimum match length for MetaCC-style filtering |
| `--metacc-min-binsize` | `150000` | Minimum bin size for MetaCC-style binning |
| `--normcc-thres` | `0.05` | NormCC threshold |
| `--bin3c-min-len` | `1000` | Minimum contig length for bin3C |
| `--bin3c-min-signal` | `5` | Minimum contact signal for bin3C |
| `--bin3c-min-mapq` | `60` | Minimum MAPQ for bin3C |
| `--bin3c-min-match` | `10` | Minimum match length for bin3C |
| `--bin3c-min-extent` | `50000` | Minimum extent for bin3C |
| `--min-completeness` | `50` | Minimum completeness for retained bins |
| `--max-contamination` | `10` | Maximum contamination for retained bins |
| `--contamination-penalty` | `5` | Penalty used in representative-bin scoring |
| `--min-input-bin-size` | `50000` | Minimum input bin size considered during refinement |
| `--max-input-bin-size` | `20000000` | Maximum input bin size considered during refinement |
| `--binning-refiner-min-size` | `524288` | Minimum size used by Binning_refiner |
| `--tmp-dir` | system temporary directory | Temporary directory root for CheckM2 working files |
| `--keep-temp` | `false` | Keep successful CheckM2 temporary directories |
| `--no-fasta` | `false` | Do not write bin3C cluster FASTA files |
| `--no-report` | `false` | Do not write bin3C cluster report |
| `--no-spades` | `true` | Mark input assembly as not produced by SPAdes |
| `--only-large` | `false` | Only write large bin3C FASTA clusters |
| `--skip-checkm2` | `false` | Skip CheckM2 during final bin refinement |
| `--skip-refinement` | `false` | Skip hybrid bin refinement |
| `--skip-consolidation` | `false` | Skip final consolidation |
| `--keep-ambiguous` | `false` | Keep ambiguous contigs in all bins |
| `--remove-ambiguous` | `false` | Remove ambiguous contigs from all bins |
| `--seed` | - | Optional random seed |

## Outputs

| Output | Description |
| --- | --- |
| `6_binning/binning/` | Binning output directory |
| `metahict/metahict_50_10_bins/` | Consolidated MAG FASTA files |
| CheckM2 reports | Completeness and contamination estimates |

## Nextflow command

```bash
nextflow run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module6 \
  --samplesheet nextflow/assets/test_data_samplesheet.csv \
  --assembly_dir "$PWD/results/test_data/test_data/2_assembly" \
  --alignment_dir "$PWD/results/test_data/test_data/3_alignment" \
  --out_root "$PWD/results/module6" \
  --report_dir "$PWD/results/module6/nextflow_reports" \
  -work-dir "$PWD/results/module6/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
export CHECKM2DB="$PWD/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"

python metahict.py binning \
  --project_path "$PWD/modules" \
  --fasta results/manual/module2/final_assembly.fasta \
  --bam results/manual/module3/sorted_map.bam \
  --output results/manual/module6 \
  --enzyme "Sau3AI,MluCI" \
  -t 80 \
  --no-spades \
  --checkm_db "$PWD/databases/checkm_db" \
  --checkm2_db "$CHECKM2DB"
```
