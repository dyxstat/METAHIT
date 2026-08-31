# Reassembly

Reassembly selects short-insert, shotgun-like intra-contig Hi-C read pairs
using an expectation-maximization (EM) model, recruits reads to each MAG,
combines them with shotgun evidence, and runs per-bin SPAdes assemblies. It can
improve continuity, but it can also introduce or reveal conflicts. METAHICT
evaluates each candidate reassembly before selecting an output.

This stage is defined for paired short-read shotgun samples. The complete
workflow skips it for long-read samples and sends the final binning MAGs to
scaffolding, annotation, and MGE–host pair analysis.

## Before you run

All upstream results must come from the same sample and compatible
assembly. Check that `6_binning/metahict/final_bins/` contains the intended
MAGs, the alignment BAM uses the original assembly, and the CheckM2 database
is available. Reassembly is CPU-, memory-, temporary-disk-, and file-intensive.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--binning-dir` | `./metahict run` | Binning output directory with consolidated bins |
| `--assembly-dir` | `./metahict run` | Original assembly directory |
| `--alignment-dir` | `./metahict run` | Hi-C alignment directory |
| `--preprocessing-dir` | `./metahict run` | Preprocessing output directory containing cleaned reads in both `sg/` and `hic/` |

## Optional parameters

Set `resources.reassembly.memory` in `metahict_configuration.yaml`; the default
is `64 GB`, and 80% becomes the total SPAdes memory budget. Set algorithm
options under `modules.reassembly` in the same file. Pass the
database path through the top-level `--checkm2-db` option.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `read_selection.em_cutoff_quantile` | `0.95` | Quantile defining the short-insert cutoff |
| `read_selection.em_top_contigs` | `100` | Longest contigs used for EM fitting |
| `read_selection.min_mapping_quality` | `30` | Minimum MAPQ for insert-size extraction |
| `read_selection.min_aligned_length` | `30` | Minimum aligned match length |
| `read_selection.exclude_duplicate_alignments` | `false` | Exclude duplicate-marked alignments |
| `read_selection.write_nonselected_hic_reads` | `false` | Write Hi-C pairs not selected by the EM cutoff |
| `read_recruitment.strict_mismatch_cutoff` | `2` | Strict mismatch cutoff for read recruitment |
| `read_recruitment.permissive_mismatch_cutoff` | `5` | Permissive mismatch cutoff for read recruitment |
| `assembly.min_contig_length` | `500` | Minimum contig length retained after reassembly |
| `assembly.spades_mode` | `careful` | SPAdes mode; use `none` to disable careful mode |
| `assembly.phred_offset` | - | SPAdes PHRED offset |
| `assembly.extra_args` | - | Additional SPAdes arguments |
| `assembly.assemble_residual_reads` | `true` | Assemble residual non-bin reads |
| `quality_control.min_completeness` | `50` | Minimum completeness used to choose a reassembly |
| `quality_control.max_contamination` | `10` | Maximum contamination used to choose a reassembly |
| `quality_control.contamination_penalty` | `5` | Penalty used to score a reassembly |
| `quality_control.run_checkm2` | `true` | Evaluate reassemblies with CheckM2 |
| `--checkm2-db` | `databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd` | CheckM2 database path passed to `./metahict run` |
| `temporary_directory` | system temporary directory | Temporary directory root |
| `keep_temporary_files` | `false` | Keep read-selection, reassembly, SPAdes, and CheckM2 intermediate files for debugging |

## Outputs

| Output | Description |
| --- | --- |
| `7_reassembly/reassembled_bins/` | Reassembled MAG FASTA files |
| `7_reassembly/reassembled_bins_quality.tsv` | Final CheckM2 quality summary for the selected MAGs |
| `7_reassembly/reassembled_bin_name_map.tsv` | Mapping from selected assembly names to simplified bin names |
| `7_reassembly/combined_contigs.fa` | Combined binned and residual contigs for downstream modules |
| `7_reassembly/residual_contigs.fa` | Assembly of residual reads not assigned to a MAG |
| `7_reassembly/read_selection_summary.json` | EM model, insert-size counts, cutoff, and read-selection diagnostics |
| `7_reassembly/quality/` | Original-bin comparison and final raw CheckM2 table |
| `7_reassembly/figures/` | Final quality and before/after comparison figures |
| `7_reassembly/run_parameters.yaml` | Effective parameters used by the reassembly run |
| `7_reassembly/intermediates/` | Large debugging files, present only when `keep_temporary_files: true` |

## How to check the result

Use `reassembled_bin_name_map.tsv` to match new names to original MAGs. Compare
`reassembled_bins_quality.tsv` with `quality/original_bins_quality.tsv`,
including any change in completeness, contamination, and sequence continuity.
Review `read_selection_summary.json` for the fitted insert-size cutoff and
number of selected Hi-C pairs. Confirm that `combined_contigs.fa` contains the
expected selected and residual sequence before MGE analysis.

## Recommended command

```bash
./metahict run \
  --entry-module reassembly \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --binning-dir results/sample_01/6_binning \
  --assembly-dir results/sample_01/2_assembly \
  --alignment-dir results/sample_01/3_alignment \
  --preprocessing-dir results/sample_01/1_preprocessing \
  --outdir results/reassembly
```

[Previous: Binning](binning.md) · [Next: Scaffolding](scaffolding.md) · [All modules](README.md)
