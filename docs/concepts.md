# Concepts for new users

This introduction is written for readers familiar with sequencing data but new to
metagenomic Hi-C.

## The two read libraries have different jobs

METAHICT expects two libraries from the same biological sample.

| Library | What the read pairs mainly tell METAHICT | Main uses |
| --- | --- | --- |
| Shotgun metagenome | Which DNA sequences are present and how deeply they were sequenced | Assembly, contig coverage, read recruitment, and reassembly |
| Metagenomic Hi-C | Which DNA fragments were physically close when the sample was cross-linked | Contig contacts, genome binning, scaffolding, and candidate MGE–host pairs |

Shotgun reads provide the primary sequence information. Hi-C read pairs are
often chimeric by design: their two ends can originate from different DNA
fragments that were ligated because they were close in the cross-linked
sample. METAHICT treats the libraries differently during analysis.

The libraries should represent the same community. With unmatched samples,
coverage and contact evidence describe different populations and cannot be
interpreted together.

## Important terms

| Term | Plain-language meaning |
| --- | --- |
| Contig | A continuous sequence assembled from overlapping reads |
| Contacts | Hi-C read pairs whose ends align to two contigs or two positions |
| Contact matrix | A table recording contact strength between contig pairs |
| Binning | Grouping contigs predicted to originate from the same genome |
| MAG | A metagenome-assembled genome: a bin reconstructed from community sequencing rather than an isolated culture |
| Reassembly | Recruiting reads for one bin and assembling them again |
| Scaffolding | Ordering and orienting contigs within a selected bin, usually leaving sequence gaps |
| MGE | A mobile genetic element, such as a viral or plasmid sequence |
| Restriction enzyme | The enzyme used to cut DNA during Hi-C library preparation; it affects where contacts can be observed |

Completeness and contamination are model-based quality estimates, not direct
measurements. Thresholds such as 50% completeness and 10% contamination are
workflow filters, not universal definitions of a good genome.

## What the ten stages do

| Stage | Question answered | Primary result |
| --- | --- | --- |
| 1. Preprocessing | Are the short reads clean enough to use? | Clean paired reads and FastQC reports; long shotgun reads bypass this stage |
| 2. Assembly | What longer sequences can be reconstructed from shotgun reads? | Contig FASTA |
| 3. Alignment | Where do the Hi-C read ends map on the contigs? | Filtered Hi-C BAM |
| 4. Coverage | How deeply is each contig represented by shotgun reads? | Contig-depth table |
| 5. Contact | Which contigs have Hi-C support, after filtering and normalization? | Raw and normalized contact matrices |
| 6. Binning | Which contigs are likely to belong to the same genome? | Consolidated MAG FASTAs |
| 7. Reassembly | Can recruited short reads improve each MAG? | Reassembled MAGs and combined contigs; skipped for long-read samples |
| 8. Scaffolding | Can Hi-C links order and orient contigs within a MAG? | Scaffolded sequence and heatmap |
| 9. Annotation | Where do the MAGs fall in the GTDB taxonomy? | GTDB-Tk summary tables |
| 10. MGE | Which contigs are candidate viruses/plasmids, and which hosts have contact support? | MGE calls and candidate MGE–host pairs |

The complete workflow carries files between these stages automatically.
Selected-module commands are useful when a stage must be repeated or when
compatible results were produced elsewhere.

## Why contacts are normalized

Raw Hi-C counts are influenced by more than biological proximity. Longer
contigs, restriction-site abundance, coverage, alignment quality, and random
or spurious ligations can all affect observed counts. The contact stage first
filters alignments and weak signals, then applies the selected normalization
method. Keep the actual restriction enzyme in each samplesheet row because
enzyme recognition sites are part of this interpretation.

Normalized contact strength measures proximity evidence. Closely related
genomes, shared or repetitive sequence, extracellular DNA, high abundance,
and experimental noise can create ambiguous links. Interpret contacts
together with bin quality, scaffold heatmaps, raw support, and independent
evidence.

## Decisions to make before the first run

1. Confirm that shotgun and Hi-C files belong to the same sample.
2. Obtain the exact restriction enzyme or enzyme list from the library
   preparation record.
3. Use short sample identifiers containing letters, numbers, underscores, or
   hyphens; do not encode experimental conclusions in the name.
4. Review any automatic CPU or memory cap warning and confirm that the
   effective resources remain suitable for memory-intensive stages.
5. Decide where large databases, temporary work files, and final results will
   be stored.
6. Use the distributed scientific defaults for the first run, then adjust
   parameters when the data or experimental design provides a reason.

Next: [Installation](installation.md) and the [Beginner quick
start](quickstart.md).
