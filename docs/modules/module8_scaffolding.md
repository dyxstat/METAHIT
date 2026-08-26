# Module 8: scaffolding

Scaffolding uses Hi-C contacts to scaffold MAG contigs with YaHS and produces contact heatmaps for inspection.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--reassembly_dir` | Nextflow selected-module run | Reassembly output directory |
| `--alignment_dir` | Nextflow selected-module run | Hi-C alignment output directory |
| `--hic_preprocessing_dir` | Nextflow selected-module run | Cleaned Hi-C reads |
| `enzyme` | Nextflow sample sheet | Restriction enzyme list |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `--fasta` | module wrapper | Input bin FASTA |
| `--bam` | module wrapper | Optional Hi-C alignment BAM |
| `--enzyme` | module wrapper | Restriction enzyme list |
| `--outdir` | module wrapper | Output directory |
| `--hic1`, `--hic2` | module wrapper | Preprocessed Hi-C paired-end reads |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--resolution` | `10000` | Heatmap bin size in bp |
| `--min-contig-len` | `5000` | Minimum contig length for scaffolding |
| `--bwa-options` | `-5SP` | BWA-MEM options |
| `--samtools-filter` | `-F 0x900` | SAM flag filter |
| `--sort-memory` | `1G` | SAMtools sort memory |
| `--metacc-min-mapq` | `30` | Minimum MAPQ for MetaCC-style filtering |
| `--metacc-min-len` | `1000` | Minimum contig length for MetaCC-style filtering |
| `--metacc-min-match` | `30` | Minimum match length for MetaCC-style filtering |
| `--metacc-min-signal` | `2` | Minimum contact signal for MetaCC-style filtering |
| `--bin3c-min-mapq` | `60` | Minimum MAPQ for bin3C-style filtering |
| `--bin3c-min-len` | `1000` | Minimum contig length for bin3C-style filtering |
| `--bin3c-min-match` | `10` | Minimum match length for bin3C-style filtering |
| `--bin3c-min-signal` | `5` | Minimum contact signal for bin3C-style filtering |
| `--yahs-resolutions` | - | Optional YaHS resolution list |
| `--yahs-min-mapq` | `10` | Minimum MAPQ for YaHS |
| `--yahs-min-contig-len` | `0` | Minimum contig length passed to YaHS |
| `--yahs-rounds` | `1` | Number of YaHS rounds |
| `--yahs-no-contig-ec` | `false` | Disable YaHS contig error correction |
| `--yahs-no-scaffold-ec` | `false` | Disable YaHS scaffold error correction |
| `--yahs-no-mem-check` | `false` | Disable YaHS runtime memory check |
| `--yahs-extra-args` | - | Extra YaHS arguments |
| `--normcc-thres` | `0.05` | NormCC threshold |
| `--heatmap-max-image` | `5000` | Maximum rendered heatmap image size |
| `--memory` | 80% of available memory | Memory limit |
| `--skip-checkm2` | `false` | Skip CheckM2 quality evaluation |
| `--checkm2_db` | `databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` | CheckM2 database path |
| `--tmp-dir` | system temporary directory | Temporary directory root |
| `--keep-temp` | `false` | Keep temporary files for debugging |

## Outputs

| Output | Description |
| --- | --- |
| `8_scaffolding/scaffolding/` | Scaffolding output directory |
| scaffold FASTA files | YaHS scaffolded MAG assemblies |
| contact heatmaps | Hi-C heatmaps for scaffold inspection |
| CheckM2 reports | Optional quality reports unless skipped |

## Nextflow command

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module8 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --reassembly_dir "$PWD/results/example_data/example_data/7_reassembly/reassembly" \
  --alignment_dir "$PWD/results/example_data/example_data/3_alignment" \
  --hic_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module8" \
  --report_dir "$PWD/results/module8/nextflow_reports" \
  -work-dir "$PWD/results/module8/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

```bash
python metahict.py scaffolding \
  -p "$PWD" \
  --fasta results/manual/module7/reassembled_bins/bin1.fa \
  --bam results/manual/module3/sorted_map.bam \
  --enzyme "Sau3AI,MluCI" \
  --outdir results/manual/module8 \
  --hic1 results/manual/module1_hic/final_*_1.fastq.gz \
  --hic2 results/manual/module1_hic/final_*_2.fastq.gz \
  -t 80 \
  --checkm2_db "$PWD/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
```
