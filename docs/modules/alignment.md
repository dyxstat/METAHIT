# Alignment

Alignment maps cleaned Hi-C read ends to the assembled contigs with BWA-MEM
and filters/sorts the alignments for contact construction. This BAM represents
Hi-C proximity evidence. The coverage stage performs a separate shotgun-read
alignment.

## Before you run

The assembly and Hi-C reads must belong to the same sample. A reusable BAM
must match the assembly's contig names and lengths exactly. Review
preprocessing QC and ensure temporary storage can hold the BAM and sort files.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--assembly-dir` | `./metahict run` | Assembly output directory containing the contig FASTA |
| `--preprocessing-dir` | `./metahict run` | Preprocessing output directory; alignment reads its `hic/` child |

## Optional parameters

Set these stage options under `modules.alignment` in `metahict_configuration.yaml`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `bwa.options` | `-5SP` | BWA-MEM options |
| `filtering.samtools_filter` | `-F 0x900` | SAM flag filter |
| `filtering.min_mapping_quality` | `30` | Minimum mapping quality |
| `filtering.min_aligned_length` | `30` | Minimum nucleotide match length |
| `filtering.min_intra_contig_distance` | `10000` | Minimum intra-contig distance used in metrics |
| `sorting.memory_per_thread` | `1G` | Memory per SAMtools sort thread |
| `sorting.temporary_directory` | system temporary directory | Temporary directory root |
| `sorting.keep_sam` | `false` | Keep intermediate SAM file |
| `metrics.enabled` | `true` | Generate alignment diagnostics |

## Outputs

| Output | Description |
| --- | --- |
| `3_alignment/sorted_map.bam` | Filtered, sorted Hi-C alignment BAM |
| alignment metrics | Alignment-derived QC metrics, including 3D ratio unless skipped |

## How to check the result

Confirm that `sorted_map.bam` and its index are non-empty, then inspect the
alignment metrics. Very low mapping or retained-pair rates may indicate an
unmatched library, poor assembly, wrong preprocessing, contamination, or
overly strict filters. Mapping quality describes confidence in placement;
contact quality also depends on the library and downstream filtering.

## Recommended command

```bash
./metahict run \
  --entry-module alignment \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --assembly-dir results/sample_01/2_assembly \
  --preprocessing-dir results/sample_01/1_preprocessing \
  --outdir results/alignment
```

[Previous: Assembly](assembly.md) · [Next: Coverage](coverage.md) · [All modules](README.md)
