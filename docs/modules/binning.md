# Binning and consolidation

Binning groups contigs into candidate metagenome-assembled genomes (MAGs).
METAHICT runs MetaCC, bin3C, and ImputeCC, constructs hybrid candidates,
evaluates them, and consolidates a non-redundant final collection.

Agreement between tools can strengthen confidence, but every MAG remains a
reconstruction. Completeness and contamination are estimates. The distributed
50% completeness and 10% contamination thresholds define the default retained
set and should be reported with the results.

## Before you run

The assembly and Hi-C BAM must use identical contig names and lengths. The
samplesheet enzyme must be correct. Verify both CheckM and CheckM2 databases,
review available memory and temporary storage, and inspect alignment and
contact QC. Poor input quality can reduce bin recovery.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--assembly-dir` | `./metahict run` | Assembly output directory containing the contig FASTA |
| `--alignment-dir` | `./metahict run` | Hi-C alignment output directory |
| `enzyme` | Samplesheet | Required restriction enzyme list |

## Optional parameters

Set algorithm options under `modules.binning` in `metahict_configuration.yaml`. Pass the
database path through the top-level `--checkm-db` option.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--checkm-db` | `databases/checkm_db` | CheckM database path passed to `./metahict run` |
| `metacc.min_contig_length` | `1000` | Minimum contig length for MetaCC |
| `metacc.min_contact_signal` | `2` | Minimum contact signal for MetaCC |
| `metacc.min_mapping_quality` | `30` | Minimum alignment MAPQ for MetaCC |
| `metacc.min_aligned_length` | `30` | Minimum aligned match length for MetaCC |
| `metacc.min_bin_size` | `150000` | Minimum MetaCC bin size |
| `metacc.normcc_discard_fraction` | `0.05` | Fraction of weakest NormCC contacts discarded |
| `metacc.marker_gene_count` | auto | Override the marker-gene count used by MetaCC |
| `bin3c.min_contig_length` | `1000` | Minimum contig length for bin3C |
| `bin3c.min_contact_signal` | `5` | Minimum contact signal for bin3C |
| `bin3c.min_mapping_quality` | `60` | Minimum alignment MAPQ for bin3C |
| `bin3c.min_aligned_length` | `10` | Minimum aligned match length for bin3C |
| `bin3c.min_cluster_extent` | `50000` | Minimum bin3C cluster extent |
| `imputecc.gene_coverage` | `0.9` | Marker-gene coverage threshold |
| `imputecc.random_walk_restart_probability` | `0.5` | Random-walk restart probability |
| `imputecc.random_walk_threshold` | `80` | Random-walk threshold |
| `imputecc.max_markers` | `8000` | Maximum marker count |
| `imputecc.intra_bin_threshold` | `50` | Intra-bin contact threshold |
| `imputecc.inter_bin_threshold` | `0` | Inter-bin contact threshold |
| `imputecc.min_bin_size` | `100000` | Minimum ImputeCC bin size |
| `imputecc.contamination_weight` | `2` | ImputeCC contamination weight |
| `imputecc.min_completeness` | `50` | Minimum completeness reported by ImputeCC |
| `imputecc.max_contamination` | `10` | Maximum contamination reported by ImputeCC |
| `imputecc.report_quality_threshold` | `10` | ImputeCC report-quality threshold |
| `refinement.min_completeness` | `50` | Minimum completeness for retained consolidated bins |
| `refinement.max_contamination` | `10` | Maximum contamination for retained consolidated bins |
| `refinement.contamination_penalty` | `5` | Penalty used in representative-bin scoring |
| `refinement.min_input_bin_size` | `50000` | Minimum input bin size considered during refinement |
| `refinement.max_input_bin_size` | `20000000` | Maximum input bin size considered during refinement |
| `refinement.binning_refiner_min_size` | `524288` | Minimum size passed to Binning_refiner |
| `refinement.run_checkm2` | `true` | Evaluate candidate bins with CheckM2 |
| `refinement.run_refinement` | `true` | Run hybrid-bin refinement |
| `refinement.run_consolidation` | `true` | Build the final non-redundant bin set |
| `refinement.ambiguous_contigs` | `best` | Handle ambiguous contigs with `best`, `keep`, or `remove` |
| `output.generate_bin3c_fasta` | `true` | Write bin3C cluster FASTA files |
| `output.generate_bin3c_report` | `true` | Write the bin3C cluster report |
| `output.input_assembly_is_spades` | `false` | Interpret SPAdes-style contig names |
| `output.only_large_bin3c_clusters` | `false` | Write only large bin3C FASTA clusters |
| `output.heatmap_max_image_size` | `5000` | Maximum heatmap image dimension |
| `random_seed` | - | Optional random seed |
| `temporary_directory` | system temporary directory | Temporary directory root |
| `keep_temporary_files` | `false` | Keep refinement and CheckM2 intermediate files for debugging |

## Outputs

| Output | Description |
| --- | --- |
| `6_binning/metahict/final_bins/` | Consolidated MAG FASTA files; this is the downstream bin input |
| `6_binning/metahict/final_bins_quality.tsv` | Completeness and contamination estimates for the final MAGs |
| `6_binning/metahict/contig_to_bin.tsv` | Final contig-to-MAG assignments |
| `6_binning/metahict/combined_final_bins.fa` | All final MAG sequences in one FASTA file |
| `6_binning/metahict/run_parameters.yaml` | Effective parameters used by the binning run |
| `6_binning/metahict/figures/` | Binning-quality and contact-map figures |
| `6_binning/metahict/intermediates/` | Debugging artifacts, present only when `keep_temporary_files: true` |

## How to check the result

Review `final_bins_quality.tsv` and the FASTAs under `final_bins/`. Inspect the
completeness and contamination distributions and check the unbinned and
ambiguous contigs. MAG selection should use these quality measures together
with the biological objective. The thresholds and other effective settings
are recorded in `run_parameters.yaml`.

## Recommended command

```bash
./metahict run \
  --entry-module binning \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --assembly-dir results/sample_01/2_assembly \
  --alignment-dir results/sample_01/3_alignment \
  --outdir results/binning
```

## CheckM2 and DIAMOND validation

CheckM2 calls DIAMOND and other companion executables from its own locked
environment. METAHICT adds that environment to `PATH`. To check
the installed DIAMOND executable independently before a long binning run, use:

```bash
conda run -p conda_envs/checkm2 diamond version
```

If this command reports that DIAMOND is missing, the `checkm2` environment is
incomplete and should be repaired before rerunning binning. A failed METAHICT
run can be repeated with the same arguments plus `--resume`; do not delete its
`nextflow_work` directory.

[Previous: Contact](contact.md) · [Next: Reassembly](reassembly.md) · [All modules](README.md)
