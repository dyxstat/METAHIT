# Assembly

Assembly reconstructs longer contigs from cleaned reads. These contigs become
the common reference for Hi-C alignment, shotgun coverage, contact matrices,
binning, and downstream contig identifiers. MEGAHIT is the default for normal
paired short-read METAHICT data; metaSPAdes is an alternative short-read
assembler.

When `long_read_type` is present in the samplesheet, METAHICT passes the single
`sg1` file directly to metaFlye. Short-read rows use the configured MEGAHIT or
metaSPAdes method.

## Before you run

Check the preprocessing FastQC reports, confirm sufficient temporary disk
space, and review the task memory. Assembly is one of the largest tasks. The
selected preprocessing directory must directly contain the cleaned pair named
`final_*_1.fastq.gz` and `final_*_2.fastq.gz`. Long-read mode instead requires
a readable single `sg1` file and a correct `long_read_type` value.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `--preprocessing-dir` | `./metahict run` | Required for short reads; assembly reads its `sg/` child |
| `sg1`, `long_read_type` | Samplesheet | Single long-read file and its metaFlye input type; `sg2` must be empty |

Valid long-read types are `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`,
`nano-raw`, `nano-corr`, and `nano-hq`.

## Optional parameters

Set `resources.assembly.memory` in `metahict_configuration.yaml`; the default is
`64 GB`. METAHICT passes 80% of that allocation to the assembler. Set the
remaining stage options under `modules.assembly` in the same file.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `assembler` | `megahit` | Short-read assembler (`megahit` or `metaspades`); long-read rows automatically use metaFlye |
| `min_contig_length` | `1000` | Minimum contig length retained |
| `temporary_directory` | system temporary directory | Temporary directory root |
| `keep_temporary_files` | `false` | Keep temporary assembler files |
| `quality_control.run_quast` | `true` | Generate a QUAST assembly report |
| `megahit.min_kmer_length` | `21` | MEGAHIT minimum k-mer size |
| `megahit.max_kmer_length` | `141` | MEGAHIT maximum k-mer size |
| `megahit.kmer_step` | `12` | MEGAHIT k-mer step size |
| `megahit.merge_level` | `20,0.95` | MEGAHIT merge-level setting |
| `megahit.extra_args` | - | Additional MEGAHIT arguments |
| `metaspades.kmer_lengths` | `21,33,55` | metaSPAdes k-mer list |
| `metaspades.extra_args` | - | Additional metaSPAdes arguments |
| `metaflye.extra_args` | - | Additional metaFlye arguments |

Only the subsection matching `assembler` is passed to the selected program;
settings for the other assemblers are ignored for that run. For a long-read
row, the samplesheet type selects metaFlye regardless of the short-read
`assembler` setting.

## Outputs

| Output | Description |
| --- | --- |
| `2_assembly/final_assembly.fasta` | Assembly FASTA used by downstream modules |
| assembler output files | Native assembler outputs |
| QUAST report | Assembly QC report unless skipped |

## How to check the result

`final_assembly.fasta` must be non-empty and have a `.fai` index. Review the
QUAST report, contig count, total assembly length, and length distribution.
High fragmentation can increase contact-matrix size and make binning harder.
Use the QC results to assess downstream suitability.

## Recommended command

```bash
./metahict run \
  --entry-module assembly \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --preprocessing-dir results/sample_01/1_preprocessing \
  --outdir results/assembly
```

The assembly is published at
`results/assembly/sample_01/2_assembly/final_assembly.fasta`.

For a long-read samplesheet, omit `--preprocessing-dir`; assembly reads `sg1`
directly:

```bash
./metahict run \
  --entry-module assembly \
  --samplesheet long_read_samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results/assembly
```

[Previous: Preprocessing](preprocessing.md) · [Next: Alignment](alignment.md) · [All modules](README.md)
