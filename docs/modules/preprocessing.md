# Preprocessing

Preprocessing cleans paired-end short-read shotgun and Hi-C reads with BBTools/BBDuk and
produces FastQC reports before and after trimming. Shotgun reads supply the
sequence and abundance information used by assembly and coverage. Hi-C reads
supply proximity links used by alignment, contact analysis, binning, and
scaffolding. For a long-read shotgun sample, the single shotgun file bypasses
this stage and only the paired Hi-C library is preprocessed.

## When to run this stage

Run preprocessing first for raw FASTQ files. The complete workflow processes
both libraries. Use a selected run to repeat QC or trimming, or to process one
library separately.

Before running, confirm that R1 and R2 are paired correctly, the files pass
`gzip -t`, and the two libraries belong to the same biological sample.

## Required parameters

| Parameter | Used by | Meaning |
| --- | --- | --- |
| `sg1`, `sg2` | Samplesheet | Shotgun paired-end FASTQ files; long-read mode uses only `sg1` |
| `hic1`, `hic2` | Samplesheet | Hi-C paired-end FASTQ files |
| `long_read_type` | Samplesheet | Empty for short reads; required for a single long-read `sg1` file |

## Optional parameters

Set these shared keys under `modules.preprocessing` in
`metahict_configuration.yaml`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `libraries.shotgun.enabled` | `true` | Process the shotgun pair |
| `libraries.shotgun.deduplicate` | `false` | Apply Clumpify duplicate removal to shotgun reads |
| `libraries.hic.enabled` | `true` | Process the Hi-C pair |
| `libraries.hic.deduplicate` | `true` | Apply Clumpify duplicate removal to Hi-C reads |
| `output_prefix` | derived from input name | Custom output-file prefix |
| `trimming.min_read_length` | `50` | Minimum read length retained after trimming |
| `trimming.quality_trim_threshold` | `10` | Quality-trimming threshold |
| `trimming.quality_trim_ends` | `r` | BBDuk quality-trimming direction |
| `trimming.force_trim_left` | `10` | Number of bases trimmed from the left end |
| `trimming.length_modulo` | `5` | BBDuk modulo trimming value |
| `trimming.adapter_trim_end` | `r` | Adapter-trimming direction |
| `trimming.adapter_kmer_length` | `23` | Adapter k-mer size |
| `trimming.min_adapter_kmer_length` | `11` | Minimum adapter k-mer size |
| `trimming.adapter_hamming_distance` | `1` | Adapter k-mer Hamming distance |
| `trimming.adapter_reference` | - | Adapter FASTA; unset uses the Conda-installed BBTools reference |
| `quality_control.run_before_trimming` | `true` | Generate FastQC reports for input reads |
| `quality_control.run_after_trimming` | `true` | Generate FastQC reports for cleaned reads |

Nextflow derives the Java heap passed to BBDuk and Clumpify from 80% of
`resources.preprocessing.memory`.

Nextflow allocates `32 GB` to each preprocessing task by default and passes 80%
of that allocation to every BBTools Java command. Thus the default Java heap is
`-Xmx25g`. Change both the task allocation and derived heap by editing
`resources.preprocessing.memory` in `metahict_configuration.yaml`, or override
it for one invocation with `-m/--memory`.

By default, shotgun deduplication is off because abundance and natural
duplicate structure may be informative. Hi-C deduplication is on to reduce
amplification-driven contact inflation. Set
`libraries.hic.deduplicate: false` when the protocol calls for retaining Hi-C
duplicates.

## Outputs

| Output | Description |
| --- | --- |
| `1_preprocessing/sg/final_*_1.fastq.gz`, `1_preprocessing/sg/final_*_2.fastq.gz` | Cleaned shotgun paired-end reads |
| `1_preprocessing/hic/final_*_1.fastq.gz`, `1_preprocessing/hic/final_*_2.fastq.gz` | Cleaned Hi-C paired-end reads |
| FastQC reports | Pre- and post-cleaning QC reports unless skipped |

The published files are directly inside `sg/` and `hic/`; there is no second
`preprocessing/` directory. Downstream stages also accept results created by
older versions with the nested layout.

Long-read samples do not produce `1_preprocessing/sg/`; their original `sg1`
file is passed directly to assembly and coverage.

## How to check the result

Confirm that both final files in each enabled library are non-empty and inspect
the before/after FastQC reports. Investigate unexpected read loss, persistent
adapter content, severe quality decline, and pair-count inconsistencies before
assembly.

## Recommended command

```bash
./metahict run \
  --entry-module preprocessing \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results/preprocessing
```

For this sample, cleaned reads are written below
`results/preprocessing/sample_01/1_preprocessing/`.

[Next: Assembly](assembly.md) · [All modules](README.md)
