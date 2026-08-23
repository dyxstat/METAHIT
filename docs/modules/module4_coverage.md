# Module 4: coverage

Coverage maps shotgun reads to assembled contigs and summarizes contig depth for binning.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- |
| `--assembly_dir` | Nextflow selected-module run | Assembly output directory containing the contig FASTA |
| `--sg_preprocessing_dir` | Nextflow selected-module run | Shotgun preprocessing output directory |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `-1` | module wrapper | Forward shotgun read file |
| `-2` | module wrapper | Reverse shotgun read file |
| `-r`, `--reference` | module wrapper | Assembly FASTA used as the mapping reference |
| `-o`, `--output` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--percent-identity` | `97` | Minimum alignment identity for BBMap |
| `--min-mapq` | `0` | Minimum mapping quality for depth counting |
| `--weight-mapq` | `0` | MAPQ weighting option for depth calculation |
| `--include-edge-bases` | `false` | Include read-length edge bases in depth and variance |
| `--max-edge-bases` | `75` | Edge-base value used when edge bases are excluded |
| `--min-contig-length` | `0` | Minimum contig length for coverage reporting |
| `--min-contig-depth` | `0` | Minimum contig depth for coverage reporting |
| `--bbmap-extra-args` | - | Additional BBMap options |
| `--memory` | unset | Java heap for BBMap, for example `64g` |
| `--tmp-dir` | system temporary directory | Temporary directory root |
| `--keep-sam` | `false` | Keep intermediate SAM file |
| `--keep-temp` | `false` | Keep temporary working and index files |

## Outputs

| Output | Description |
| --- | --- |
| `4_coverage/coverage/` | Coverage output directory |
| contig-depth table | Coverage table used by downstream binning |

## Nextflow command

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module4 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --sg_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/sg" \
  --out_root "$PWD/results/module4" \
  --report_dir "$PWD/results/module4/nextflow_reports" \
  -work-dir "$PWD/results/module4/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py coverage \
  -p "$PWD" \
  -r results/manual/module2/final_assembly.fasta \
  -1 results/manual/module1_sg/final_*_1.fastq.gz \
  -2 results/manual/module1_sg/final_*_2.fastq.gz \
  -o results/manual/module4 \
  -t 80
```
