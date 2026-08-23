# Module 9: annotation

Annotation classifies MAGs with GTDB-Tk using GTDB release 220 by default.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--reassembly_dir` | Nextflow selected-module run | Reassembly output directory containing reassembled bins |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `--bin` | module wrapper | Directory containing input bins |
| `--outdir` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--extension` | `fa` | Genome file extension processed by GTDB-Tk |
| `--prefix` | `gtdbtk` | Prefix for GTDB-Tk output files |
| `--pplacer-cpus` | same as thread count | CPU count for pplacer |
| `--gtdbtk_db`, `--gtdbtk-db` | `databases/gtdbtk_db/release220` | GTDB-Tk database path |
| `--skip-ani-screen` | `true` | Skip GTDB-Tk ANI screening |
| `--no-skip-ani-screen` | disabled by default | Enable GTDB-Tk ANI/Mash screening; requires `--mash-db` |
| `--mash-db` | - | Mash database path when ANI screening is enabled |
| `--min-perc-aa` | `10` | Minimum percentage of amino acids in the MSA |
| `--min-af` | `0.5` | Minimum alignment fraction for species assignment |
| `--full-tree` | `false` | Use GTDB-Tk full bacterial tree |
| `--scratch-dir` | - | Scratch directory for pplacer disk-backed memory reduction |
| `--tmp-dir` | system temporary directory | Temporary directory root |
| `--force` | `false` | Continue if GTDB-Tk errors on one genome |
| `--keep-intermediates` | `false` | Keep GTDB-Tk intermediate files |
| `--debug` | `false` | Keep GTDB-Tk debug intermediates |
| `--write-single-copy-genes` | `false` | Output unaligned single-copy marker genes |
| `--gtdbtk-extra-args` | - | Extra arguments passed to `gtdbtk classify_wf` |

## Outputs

| Output | Description |
| --- | --- |
| `9_annotation/annotation/` | Annotation output directory |
| `gtdbtk.*.summary.tsv` | GTDB-Tk summary tables |
| GTDB-Tk native outputs | Classification workflow outputs |

## Nextflow command

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module9 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --reassembly_dir "$PWD/results/example_data/example_data/7_reassembly/reassembly" \
  --out_root "$PWD/results/module9" \
  --report_dir "$PWD/results/module9/nextflow_reports" \
  -work-dir "$PWD/results/module9/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py annotation \
  -p "$PWD" \
  --bin results/manual/module7/reassembled_bins \
  --outdir results/manual/module9 \
  -t 80 \
  --gtdbtk_db "$PWD/databases/gtdbtk_db/release220"
```
