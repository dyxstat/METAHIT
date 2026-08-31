# Mobile genetic element analysis

The MGE module performs one complete analysis. It identifies viral and plasmid
contigs with geNomad, finds circular contigs with geNomad and ccfind, and uses
Hi-C contacts to report candidate MGE–host pairs.

The input FASTA can come from short-read, long-read, hybrid, or reassembled
metagenomic data. Host-genome FASTAs must contain contigs from the same input
FASTA, and the Hi-C evidence must have been generated against that same
sequence set.

All MGE tasks use the single `resources.mge` setting, which defaults to 16
threads and 32 GB. Contact normalization is serial and therefore uses one CPU.

## Before you run

Prepare:

- A metagenome FASTA containing both MGE candidates and host contigs.
- A directory containing one `.fa`, `.fasta`, or `.fna` file per host MAG.
- Cleaned Hi-C reads, an existing MGE alignment, or an existing MGE contact
  directory generated from the same FASTA.
- A geNomad database.
- The actual restriction enzyme in the samplesheet when METAHICT must construct
  the alignment and contact matrix.

FASTA identifiers must be unique. Every contig in `--host-dir` must match one
identifier in `--fasta`. A host FASTA may contain a bare contig identifier such
as `contig_12` when the input FASTA contains either `contig_12` or
`host_name|contig_12`.

## Required parameters

| Parameter | Meaning |
| --- | --- |
| `--samplesheet` | CSV providing the sample name and, when needed, Hi-C reads and restriction enzyme |
| `--fasta` | Metagenome FASTA containing the MGE and host contigs |
| `--host-dir` | Directory containing the host MAG FASTAs |
| `--preprocessing-dir` | Preprocessing output whose `hic/` child is used to generate MGE-specific alignment and contact data |

For a multi-sample sheet, include `{sample}` in `--fasta`, `--host-dir`, and
`--preprocessing-dir` so that each row resolves to its own files.

Instead of `--preprocessing-dir`, use `--mge-alignment-dir` to reuse an
MGE-specific alignment or `--mge-contact-dir` to reuse raw and normalized
contact matrices. Reused files must have been generated from the same FASTA.

## Configuration parameters

Set algorithm options under `modules.mge` in
`metahict_configuration.yaml`. Supply a non-default database path with
`--genomad-db`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `keep_intermediates` | `false` | Publish generated Hi-C alignment and contact directories under `10_MGE/intermediates/` |
| `--genomad-db` | `databases/genomad_db/genomad_db` | geNomad database path |
| `genomad.splits` | `8` | Number of geNomad database splits |
| `genomad.sensitivity` | `4.2` | geNomad MMseqs2 sensitivity |
| `genomad.cleanup` | `true` | Delete geNomad intermediate files |
| `genomad.restart` | `false` | Overwrite existing geNomad intermediate files |
| `genomad.preset` | `default` | geNomad filtering preset |
| `genomad.min_score` | `0.7` | Minimum geNomad score with the default preset |
| `genomad.max_false_discovery_rate` | `0.1` | Maximum geNomad false-discovery rate with the default preset |
| `genomad.extra_args` | empty | Additional geNomad arguments |
| `ccfind.terminal_fragment_size` | `500` | Length examined at each contig end |
| `ccfind.min_percent_identity` | `94` | Minimum terminal-alignment identity |
| `ccfind.min_aligned_length` | `50` | Minimum terminal-alignment length |
| `pairs.filter_method` | `zscore` | Pair filter: `zscore`, `fixed`, `percentage`, or `raw-support-only` |
| `pairs.zscore_threshold` | `0.5` | Minimum Z-score used by `zscore` filtering |
| `pairs.fixed_contact_threshold` | `0` | Minimum normalized contact used by `fixed` filtering |
| `pairs.top_percent` | `50` | Percentage retained by `percentage` filtering |
| `pairs.min_raw_contacts` | `2` | Minimum raw Hi-C read-pair support |
| `pairs.min_contact_strength` | `0` | Minimum normalized contact strength considered |
| `temporary_directory` | system temporary directory | Temporary directory root |

When the MGE module generates its own Hi-C evidence, configure the internal
steps under `modules.mge.alignment` and `modules.mge.contact`.

### Generated Hi-C alignment settings

| Parameter | Default | Meaning |
| --- | --- | --- |
| `alignment.bwa.options` | `-5SP` | BWA-MEM options |
| `alignment.filtering.samtools_filter` | `-F 0x900` | SAM flag filter |
| `alignment.filtering.min_mapping_quality` | `30` | Minimum mapping quality |
| `alignment.filtering.min_intra_contig_distance` | `10000` | Minimum intra-contig distance used in metrics |
| `alignment.filtering.min_aligned_length` | `30` | Minimum aligned match length |
| `alignment.sorting.memory_per_thread` | `1G` | Memory per SAMtools sort thread |
| `alignment.sorting.temporary_directory` | system temporary directory | Temporary directory root |
| `alignment.sorting.keep_sam` | `false` | Keep intermediate SAM output |
| `alignment.metrics.enabled` | `true` | Generate alignment diagnostics |

### Generated contact settings

| Parameter | Default | Meaning |
| --- | --- | --- |
| `contact.method` | `normcc` | Contact-normalization method |
| `contact.raw_contacts.min_contact_signal` | `1` | Minimum retained contact signal |
| `contact.raw_contacts.min_contig_length` | `1000` | Minimum contig length |
| `contact.raw_contacts.min_mapping_quality` | `30` | Minimum alignment mapping quality |
| `contact.raw_contacts.min_aligned_length` | `30` | Minimum aligned match length |
| `contact.denoising.spurious_contact_percent` | `5` | Lowest normalized-contact percentage removed |
| `contact.normcc.epsilon` | `1` | NormCC numerical stabilizer |
| `contact.hiczin.coverage_file` | unset | Matching coverage table required for HiCzin |
| `contact.hiczin.epsilon` | `1` | HiCzin numerical stabilizer |
| `contact.bin3c.epsilon` | `1` | bin3C numerical stabilizer |
| `contact.bin3c.max_iterations` | `1000` | Maximum bin3C balancing iterations |
| `contact.bin3c.convergence_tolerance` | `1e-6` | bin3C convergence tolerance |
| `contact.metator.coverage_file` | unset | Matching coverage table required for MetaTOR |
| `contact.metator.epsilon` | `1` | MetaTOR numerical stabilizer |

Only the selected contact-method subsection is used. These alignment and
contact settings are not used when `--mge-contact-dir` is supplied.

## Outputs

Final products are published directly under `10_MGE/`. Generated alignment and
contact files are supporting intermediates and are omitted by default.

| Output | Description |
| --- | --- |
| `10_MGE/input_assembly/contig_to_host.tsv` | Validated contig-to-host membership |
| `10_MGE/sequence_topology.tsv` | Per-contig MGE and circularity evidence |
| `10_MGE/candidate_mge_host_pairs_*_filtered.tsv` | Filtered candidate MGE–host pairs |
| `10_MGE/mge_reports/MGE_summary.txt` | MGE, circularity, and pair-analysis summary |
| `10_MGE/genomad_output/` | Viral and plasmid calls and taxonomy |
| `10_MGE/ccfind_output/` | Circular-contig calls |
| `10_MGE/intermediates/alignment/` | Generated MGE-specific Hi-C alignment when `keep_intermediates: true` |
| `10_MGE/intermediates/contact/` | Generated raw and normalized contact matrices when `keep_intermediates: true` |

Review the geNomad scores and circularity evidence before interpreting the
MGE–host pairs. Pair calls are candidates for follow-up, not proof that an MGE
is biologically linked to a host.

## Recommended command

```bash
./metahict run \
  --entry-module mge \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --fasta inputs/metagenome.fasta \
  --host-dir inputs/host_genomes \
  --preprocessing-dir results/preprocessing/sample_01/1_preprocessing \
  --outdir results/mge \
  --threads 16
```

## Run with METAHICT reassembly outputs

```bash
./metahict run \
  --entry-module mge \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --fasta results/reassembly/sample_01/7_reassembly/combined_contigs.fa \
  --host-dir results/reassembly/sample_01/7_reassembly/reassembled_bins \
  --preprocessing-dir results/preprocessing/sample_01/1_preprocessing \
  --outdir results/mge \
  --threads 16
```

[Previous: Annotation](annotation.md) · [All modules](README.md)
