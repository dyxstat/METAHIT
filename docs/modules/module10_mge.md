# Module 10: MGE analysis

The MGE module identifies viral and plasmid contigs with geNomad, reports candidate MGE-MAG proximity associations from Hi-C contacts, and records sequence-topology annotations from geNomad and ccfind.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- |
| `--reassembly_dir` | Nextflow selected-module run | Reassembly output directory containing `combined_contigs.fa` |
| `--hic_preprocessing_dir` | Nextflow selected-module run | Cleaned Hi-C reads if MGE-specific alignment/contact steps are needed |
| `--mge_alignment_dir` | Nextflow selected-module run | Optional pre-existing MGE alignment output |
| `--mge_contact_dir` | Nextflow selected-module run | Optional pre-existing MGE contact output |
| `-p`, `--project_path` | module wrapper | METAHICT project path |
| `--combined` | module wrapper | Combined contigs FASTA |
| `--contact` | module wrapper | Normalized contact matrix |
| `--raw-contact` | module wrapper | Raw contact matrix used for raw Hi-C support |
| `--outdir` | module wrapper | Output directory |
| `-t`, `--threads` | module wrapper | Number of CPU threads; default is `80` |

## Optional parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--genomad_db`, `--genomad-db` | `databases/genomad_db/genomad_db` | geNomad database path |
| `--genomad-splits` | `8` | Number of geNomad splits |
| `--genomad-sensitivity` | `4.2` | geNomad sensitivity |
| `--genomad-cleanup` | `true` | Delete geNomad intermediate files |
| `--no-genomad-cleanup` | disabled | Keep geNomad intermediate files |
| `--genomad-restart` | `false` | Overwrite existing geNomad intermediate files |
| `--genomad-preset` | `default` | geNomad preset |
| `--genomad-min-score` | `0.7` | Minimum geNomad score |
| `--genomad-max-fdr` | `0.1` | Maximum geNomad FDR |
| `--genomad-extra-args` | - | Extra geNomad arguments |
| `--association-filter` | `zscore` | Association filtering mode |
| `--zscore-threshold` | `0.5` | Z-score threshold for default association filtering |
| `--fixed-contact-threshold` | `0` | Fixed contact threshold when using fixed filtering |
| `--top-percent` | `50` | Percentage retained when using percentage filtering |
| `--min-raw-contacts` | `2` | Minimum raw Hi-C read-pair support |
| `--ccfind-terminal-fragment-size` | `500` | ccfind terminal-fragment size |
| `--ccfind-min-identity` | `94` | ccfind minimum percent identity |
| `--ccfind-min-aligned-length` | `50` | ccfind minimum aligned length |
| `--min-contact-strength` | `0` | Minimum normalized contact strength |
| `--tmp-dir` | system temporary directory | Temporary directory root |

## Outputs

| Output | Description |
| --- | --- |
| `10_MGE/mge/candidate_mge_mag_associations_zscore_filtered.tsv` | Filtered candidate MGE-MAG proximity associations |
| `10_MGE/mge/sequence_topology.tsv` | Per-contig MGE and circularity annotations |
| geNomad outputs | Viral and plasmid calls and taxonomy |
| ccfind outputs | Circularity calls from ccfind |


## Nextflow command

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module10 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --reassembly_dir "$PWD/results/example_data/example_data/7_reassembly/reassembly" \
  --hic_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module10" \
  --report_dir "$PWD/results/module10/nextflow_reports" \
  -work-dir "$PWD/results/module10/nextflow_work" \
  -ansi-log false
```

If MGE alignment and MGE contact results already exist, reuse them:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module10 \
  --samplesheet nextflow/assets/hg_mge_samplesheet.csv \
  --reassembly_dir "$PWD/hg/7_reassembly/results/reassembly_sg_hic" \
  --mge_contact_dir "$PWD/hg/10_MGE/results/mge_contact" \
  --out_root "$PWD/hg/10_MGE/rerun" \
  --report_dir "$PWD/hg/10_MGE/rerun/nextflow_reports" \
  -work-dir "$PWD/hg/10_MGE/rerun/nextflow_work" \
  -ansi-log false
```

## Manual wrapper

First align Hi-C reads to the combined contigs produced by reassembly:

```bash
python metahict.py alignment \
  -p "$PWD" \
  -r results/manual/module7/combined_contigs.fa \
  -1 results/manual/module1_hic/final_*_1.fastq.gz \
  -2 results/manual/module1_hic/final_*_2.fastq.gz \
  -o results/manual/module10_mge_alignment \
  -t 80
```

Then build the MGE-specific contact matrices:

```bash
python metahict.py contact normcc \
  -p "$PWD" \
  --bam results/manual/module10_mge_alignment/sorted_map.bam \
  --fasta results/manual/module7/combined_contigs.fa \
  --out results/manual/module10_mge_contact \
  --enzyme "Sau3AI,MluCI"
```

Then run MGE detection and candidate MGE-MAG association analysis:

```bash
python metahict.py mge \
  -p "$PWD" \
  --combined results/manual/module7/combined_contigs.fa \
  --contact results/manual/module10_mge_contact/denoised_contact_matrix_normcc.npz \
  --raw-contact results/manual/module10_mge_contact/Raw_contact_matrix.npz \
  --outdir results/manual/module10 \
  -t 80 \
  --genomad_db "$PWD/databases/genomad_db/genomad_db"
```
