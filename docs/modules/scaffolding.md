# Scaffolding

Scaffolding uses Hi-C contacts to scaffold one selected MAG with
YaHS and produces a contact heatmap for inspection. By default, the module
aligns cleaned Hi-C reads to that bin. It can also accept a BAM whose reference
sequences and lengths match the selected bin.

Scaffolding orders and orients contigs and may insert gaps. The result remains
a draft genome whose joins and residual gaps require review. Run it after
choosing a MAG for further structural analysis. The default complete workflow
does not run scaffolding automatically.

## Before you run

Confirm the exact one-bin FASTA, the actual restriction enzyme, and enough
Hi-C support for that bin. If supplying a BAM, its reference dictionary must
match every retained FASTA sequence and length. A whole-assembly BAM is
compatible only when its reference dictionary matches the selected bin.

At least two contigs must pass `input_filter.min_contig_length` before they can
be joined. A MAG that does not meet this condition is preserved as
`unscaffolded_bin.fa`, recorded as skipped in `scaffolding_status.tsv`, and
does not cause the workflow to fail.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--scaffolding-bin` | `./metahict run` | FASTA file for the one bin to scaffold |
| `--preprocessing-dir` | `./metahict run` | Preprocessing output directory; scaffolding reads its `hic/` child |
| `enzyme` | Samplesheet | Required restriction enzyme list; there is no scientific default |

`scaffolding.py` also requires the project path, bin FASTA, Hi-C read pair,
output directory, enzyme list, and thread count. Nextflow supplies all of these:
`--scaffolding-bin` supplies the FASTA, the cleaned-read directory supplies the
Hi-C pair, the matching samplesheet row supplies the enzyme, `--outdir`
supplies the publication root, and `resources.scaffolding.threads` controls the
normal process CPU allocation. An optional `-t/--threads` value overrides it
for that run.
Missing enzyme values fail during `./metahict` preflight instead of being
silently replaced.

## Optional parameters

Set algorithm options under `modules.scaffolding` in
`metahict_configuration.yaml`. Pass the database path through the top-level
`--checkm2-db` option.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--scaffolding-bam` | unset | Optional BAM aligned to the selected bin; when unset, METAHICT performs the initial Hi-C alignment |
| `input_filter.min_contig_length` | `5000` | Minimum bin-contig length retained for scaffolding |
| `alignment.bwa_options` | `-5SP` | BWA-MEM options |
| `alignment.samtools_filter` | `-F 0x900` | SAM flag filter |
| `alignment.sort_memory_per_thread` | `1G` | SAMtools sort memory per thread |
| `contacts.min_mapping_quality` | `30` | Minimum MAPQ for contact construction |
| `contacts.min_contig_length` | `1000` | Minimum contig length for contact construction |
| `contacts.min_aligned_length` | `30` | Minimum aligned match length |
| `contacts.min_contact_signal` | `2` | Minimum retained contact signal |
| `contacts.normcc_discard_fraction` | `0.05` | Fraction of weakest NormCC contacts discarded |
| `yahs.resolutions` | - | Optional YaHS resolution list |
| `yahs.min_mapping_quality` | `10` | Minimum MAPQ for YaHS |
| `yahs.min_contig_length` | `0` | Minimum contig length passed to YaHS |
| `yahs.rounds` | `1` | Number of YaHS rounds |
| `yahs.contig_error_correction` | `true` | Enable YaHS contig error correction |
| `yahs.scaffold_error_correction` | `true` | Enable YaHS scaffold error correction |
| `yahs.memory_check` | `true` | Enable the YaHS runtime memory check |
| `yahs.extra_args` | - | Additional YaHS arguments |
| `heatmap.segment_resolution` | `10000` | Heatmap bin size in bp |
| `heatmap.max_image_size` | `5000` | Maximum rendered heatmap image dimension |
| `quality_control.run_checkm2` | `true` | Evaluate scaffolded sequences with CheckM2 |
| `--scaffolding-skip-checkm2` | unset | One-run command-line override that skips CheckM2 and its database requirement |
| `--checkm2-db` | `databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` | CheckM2 database path passed to `./metahict run` |
| `temporary_directory` | system temporary directory | Temporary directory root |
| `keep_temporary_files` | `false` | Retain alignment, MetaCC, and YaHS working files under `intermediates/` |

All algorithm options in the table other than the identified
top-level options are YAML keys. For example:

```yaml
modules:
  scaffolding:
    input_filter:
      min_contig_length: 3000
    yahs:
      min_mapping_quality: 20
    heatmap:
      segment_resolution: 20000
```

Scaffolding memory is configured as a Nextflow process resource (`32 GB` by
default). YaHS has no total-memory command-line setting.

The scaffolding contact matrix is MetaCC-only. bin3C filtering and clustering
remain in the binning module, where their outputs are used; scaffolding no
longer constructs an unused second bin3C matrix.

## Outputs

| Output | Description |
| --- | --- |
| `8_scaffolding/<BIN_ID>/scaffolding_status.tsv` | Completed or skipped status and input-contig summary |
| `8_scaffolding/<BIN_ID>/scaffolded_bin.fa` | Final scaffolded MAG |
| `8_scaffolding/<BIN_ID>/scaffolded_bin.agp` | Final scaffold structure reported by YaHS |
| `8_scaffolding/<BIN_ID>/scaffolding_metrics.txt` | Before-and-after assembly, quality, and Hi-C metrics |
| `8_scaffolding/<BIN_ID>/figures/` | Hi-C contact heatmap in PNG and PDF formats |
| `8_scaffolding/<BIN_ID>/quality/` | Original and scaffolded CheckM2 reports, unless CheckM2 is skipped |
| `8_scaffolding/<BIN_ID>/scaffolding.log` | Scaffolding program log |
| `8_scaffolding/<BIN_ID>/intermediates/` | Optional working files, present only when `keep_temporary_files: true` |
| `8_scaffolding/<BIN_ID>/unscaffolded_bin.fa` | Original MAG, written instead of scaffold outputs when fewer than two contigs pass the length filter |

The default result for one bin is therefore:

```text
8_scaffolding/<BIN_ID>/
├── scaffolding_status.tsv
├── scaffolded_bin.fa
├── scaffolded_bin.agp
├── scaffolding_metrics.txt
├── scaffolding.log
├── figures/
└── quality/
```

`intermediates/` is omitted after a successful default run. If retained for
debugging, it contains the filtered input, alignments, contact matrices,
scaffold mapping, and raw YaHS files under one directory.

For a skipped MAG, the directory contains `scaffolding_status.tsv`,
`unscaffolded_bin.fa`, and `scaffolding.log` instead of the scaffold, AGP,
metrics, and heatmap files.

## How to check the result

Inspect `scaffolded_bin.fa`, `scaffolding_metrics.txt`, and the
heatmap together.
A useful heatmap normally shows strong local structure along scaffold
diagonals; abrupt blocks or off-diagonal signals can indicate joins needing
review. Compare CheckM2 estimates and sequence content before and after
scaffolding.

## Recommended command

```bash
./metahict run \
  --entry-module scaffolding \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --scaffolding-bin results/reassembly/sample_01/7_reassembly/reassembled_bins/BIN_NAME.fa \
  --preprocessing-dir results/preprocessing/sample_01/1_preprocessing \
  --outdir results/scaffolding \
  --threads 8
```

Replace `BIN_NAME.fa` with the bin that you want to scaffold. For a
multi-sample samplesheet, use `{sample}` in `--scaffolding-bin` when each
sample has a correspondingly structured bin path.

Confirm that the samplesheet has the actual Hi-C restriction enzyme before
running:

```bash
head -n 2 samplesheet.csv
```

For example, `enzyme` may contain `Sau3AI` or `Sau3AI,MluCI`, but it
must match the Hi-C library preparation rather than the example value.

### Optional existing BAM

To avoid the initial alignment, add a BAM aligned to the same selected bin:

```bash
./metahict run \
  --entry-module scaffolding \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --scaffolding-bin results/reassembly/sample_01/7_reassembly/reassembled_bins/BIN_NAME.fa \
  --scaffolding-bam inputs/BIN_NAME.hic.bam \
  --preprocessing-dir results/preprocessing/sample_01/1_preprocessing \
  --outdir results/scaffolding \
  --threads 8
```

METAHICT checks that every sequence retained from the bin FASTA occurs in the
BAM reference dictionary with the same length, and then name-sorts the BAM.
A whole-assembly BAM can be reused when its reference sequences match the
selected bin. Cleaned Hi-C reads remain required because scaffolding realigns
them to the new scaffold for the final MetaCC matrix and heatmap.

The bundled example test invokes this standalone entry for every MAG it
recovers so that scaffolding remains covered by functional testing. Routine
complete runs omit it. Results from an explicitly selected MAG are published
below `8_scaffolding/<BIN_ID>/`.

[Previous: Reassembly](reassembly.md) · [Next: Annotation](annotation.md) · [All modules](README.md)
