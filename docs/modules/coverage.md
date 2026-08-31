# Coverage

Coverage maps paired short reads or one long-read shotgun file to the assembly and summarizes depth for
each contig. Coverage is complementary to Hi-C contact strength: it helps
describe abundance and supplies depth information used by downstream methods.

## Before you run

Use the same assembly and biological sample used for Hi-C analysis. Contig
identifiers in `coverage.txt` must remain consistent with the assembly. Check
that the shotgun preprocessing directory contains a non-empty cleaned pair.
Long-read mode maps the samplesheet `sg1` file directly.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--assembly-dir` | `./metahict run` | Assembly output directory containing the contig FASTA |
| `--preprocessing-dir` | `./metahict run` | Required for short reads; coverage reads its `sg/` child |
| `sg1`, `long_read_type` | Samplesheet | Single long-read input and type; `sg2` must be empty |

## Optional parameters

Set `resources.coverage.memory` in `metahict_configuration.yaml`; the default is
`32 GB`. METAHICT passes 80% of that allocation to BBMap as Java heap. Set the
remaining stage options under `modules.coverage` in the same file.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `mapping.min_percent_identity` | `97` | Minimum short-read alignment identity included in the depth summary |
| `mapping.long_read_min_percent_identity` | `0` | Minimum long-read alignment identity included in the depth summary; increase only after validation for the platform and data quality |
| `mapping.bbmap_extra_args` | - | Additional BBMap options |
| `depth.min_mapping_quality` | `0` | Minimum mapping quality for depth counting |
| `depth.mapping_quality_weight` | `0` | MAPQ weighting option for depth calculation |
| `depth.include_edge_bases` | `false` | Include read-length edge bases in depth and variance |
| `depth.max_excluded_edge_bases` | `75` | Edge-base value used when edge bases are excluded |
| `depth.min_contig_length` | `0` | Minimum contig length for coverage reporting |
| `depth.min_contig_depth` | `0` | Minimum contig depth for coverage reporting |
| `temporary_directory` | system temporary directory | Temporary directory root |
| `keep_sam` | `false` | Keep intermediate SAM file |
| `keep_temporary_files` | `false` | Keep temporary working and index files |

## Outputs

| Output | Description |
| --- | --- |
| `4_coverage/coverage.txt` | Contig-depth table used by downstream analyses |
| `4_coverage/pair.txt` | Paired-read coverage summary; short-read mode only |

## How to check the result

Confirm that `coverage.txt` contains a header and contig rows, and compare its
contig identifiers with `final_assembly.fasta`. Unexpectedly low mapped depth
can reflect an unmatched sample, stringent identity setting, contamination,
or an assembly dominated by unsupported contigs.

## Recommended command

```bash
./metahict run \
  --entry-module coverage \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --assembly-dir results/sample_01/2_assembly \
  --preprocessing-dir results/sample_01/1_preprocessing \
  --outdir results/coverage
```

For a long-read samplesheet, retain `--assembly-dir` but omit
`--preprocessing-dir`; coverage maps the samplesheet `sg1` file directly.

[Previous: Alignment](alignment.md) · [Next: Contact](contact.md) · [All modules](README.md)
