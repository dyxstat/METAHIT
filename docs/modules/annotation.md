# Annotation

Annotation classifies MAGs with GTDB-Tk using GTDB release 220 by default.
The distributed starting request is 8 threads and 64 GB RAM with split-tree
classification. Some datasets require more memory, and full-tree mode requires
substantially larger site-specific resources.

This stage assigns GTDB taxonomy. Gene and pathway annotation are outside its
scope.

In a complete run, paired short-read samples use reassembled MAGs and long-read
samples use the final stage 6 MAGs.

## Before you run

The selected directory must directly contain genome FASTAs with the extension
configured by `genome_extension`. Review MAG completeness and contamination,
verify the GTDB-Tk database, and compare the 64 GB starting request with the
available host memory.
METAHICT automatically caps the request to detected local memory, but a
memory-intensive classification may still require a larger allocation and
should be tested on representative data.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--samplesheet` | `./metahict run` | CSV providing the sample name |
| `--mag-dir` | `./metahict run` | Directory containing the `.fa`, `.fasta`, or `.fna` MAGs to classify |

For a multi-sample sheet, include `{sample}` in `--mag-dir`; METAHICT replaces
it with each samplesheet value. A one-sample sheet can use a direct path.
Pass the directory containing the FASTAs through `--mag-dir`.

## Optional parameters

Set algorithm options under `modules.annotation` in `metahict_configuration.yaml`. Pass
the database path through the top-level `--gtdbtk-db` option.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `genome_extension` | `fa` | Genome file extension processed by GTDB-Tk; the selected MAG directory must contain this extension |
| `output_prefix` | `gtdbtk` | Prefix for GTDB-Tk output files |
| `pplacer_threads` | same as task threads | CPU count for pplacer; cannot exceed the resolved annotation task threads after a command-line override |
| `--gtdbtk-db` | `databases/gtdbtk_db/release220` | GTDB-Tk database path passed to `./metahict run` |
| `ani_screen.enabled` | `false` | Enable GTDB-Tk ANI/Mash screening |
| `ani_screen.mash_database` | - | Mash database required when ANI screening is enabled |
| `classification.min_amino_acid_percent` | `10` | Minimum percentage of amino acids in the MSA |
| `classification.min_alignment_fraction` | `0.5` | Minimum alignment fraction for species assignment |
| `classification.use_full_tree` | `false` | Use the GTDB-Tk full bacterial tree |
| `scratch_directory` | - | Scratch directory for pplacer disk-backed memory reduction |
| `temporary_directory` | system temporary directory | Temporary directory root |
| `continue_on_genome_error` | `false` | Continue if GTDB-Tk errors on one genome |
| `keep_intermediates` | `false` | Keep GTDB-Tk intermediate files |
| `debug` | `false` | Keep GTDB-Tk debug intermediates |
| `write_single_copy_genes` | `false` | Output unaligned single-copy marker genes |
| `extra_args` | - | Additional arguments passed to `gtdbtk classify_wf` |

## Outputs

| Output | Description |
| --- | --- |
| `9_annotation/` | Annotation output directory |
| `gtdbtk.*.summary.tsv` | GTDB-Tk summary tables |
| GTDB-Tk native outputs | Classification workflow outputs |

## How to check the result

Review the bacterial and archaeal `gtdbtk.*.summary.tsv` files. Check
warnings and failed-genome records before treating a missing classification as
absence. Record the GTDB release because taxonomy can change between database
releases.

## Recommended command

```bash
./metahict run \
  --entry-module annotation \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --mag-dir results/reassembly/sample_01/7_reassembly/reassembled_bins \
  --outdir results/annotation \
  --threads 8
```

For a multi-sample sheet:

```bash
./metahict run \
  --entry-module annotation \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --mag-dir results/reassembly/{sample}/7_reassembly/reassembled_bins \
  --outdir results/annotation \
  --threads 8
```

[Previous: Scaffolding](scaffolding.md) · [Next: MGE](mge.md) · [All modules](README.md)
