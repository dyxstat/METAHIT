# Configuration

METAHICT separates sample facts from analysis choices:

- `samplesheet.csv` contains sample names, read-file locations, shotgun read
  type, and Hi-C restriction enzymes;
- `metahict_configuration.yaml` contains resource allocations and scientific
  parameters shared by the workflow.

This makes the samplesheet short and allows one reviewed configuration to be
used consistently across samples.

## Create an editable configuration

```bash
./metahict config
```

This copies the maintained template from
`nextflow/assets/metahict_configuration.yaml` to
`metahict_configuration.yaml`. An existing output is preserved unless
`--force` is supplied.

Pass the file to every run:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results
```

METAHICT records the configuration's absolute path and SHA-256 checksum in
the immutable run record.

## How the YAML file is organized

There are two top-level maps:

```yaml
resources:
  preprocessing: {threads: 8, memory: "32 GB"}
  assembly: {threads: 16, memory: "64 GB"}

modules:
  preprocessing:
    libraries:
      shotgun:
        enabled: true
        deduplicate: false
      hic:
        enabled: true
        deduplicate: true
  assembly:
    assembler: megahit
    min_contig_length: 1000
```

- `resources.<stage>` controls the CPUs and memory reserved for one task of
  that stage.
- `modules.<stage>` controls its scientific behavior.

Nested names identify the method and purpose of each setting. For example,
alignment filtering uses
`modules.alignment.filtering.min_mapping_quality`.

## YAML syntax

- Indentation is meaningful. Use spaces, not tabs.
- Write booleans as unquoted `true` or `false`.
- Write numbers without quotes unless the template already treats the value
  as text, such as `"20,0.95"`.
- Keep memory values quoted, for example `"64 GB"`.
- `null` means no explicit value; it is not the text `"null"`.
- Preserve the spelling and nesting of keys from the template.
- Keep a copy of the configuration used for a published analysis.

METAHICT rejects unknown and obsolete keys during preflight. This prevents a
misspelling from being silently ignored. A partial configuration may omit
values it does not change; workflow defaults fill those values. For clarity
and provenance, the generated complete template is recommended.

## Default task resources

| Stage | Threads | Memory | Why it is set this way |
| --- | ---: | ---: | --- |
| Preprocessing | 8 | 32 GB | Clumpify can retain large read sets in memory |
| Assembly | 16 | 64 GB | Short-read assemblers are CPU- and memory-intensive |
| Alignment | 16 | 32 GB | BWA and SAMtools use parallel threads |
| Coverage | 16 | 32 GB | BBMap uses a substantial Java heap |
| Contact | 1 | 32 GB | Current normalization is serial; matrix memory scales with contig count |
| Binning | 16 | 64 GB | Multiple binners, refinement, and CheckM2 share the task |
| Reassembly | 16 | 64 GB | Parallel per-bin SPAdes jobs share a total budget |
| Scaffolding | 8 | 32 GB | Each task scaffolds one selected bin |
| Annotation | 8 | 64 GB | GTDB-Tk classification is memory-intensive and may require an increase |
| MGE | 16 | 32 GB | One allocation covers MGE detection, circular-contig discovery, and MGE–host pairing |

These allocations are starting values. Deep read sets and large or fragmented
assemblies can need more memory and disk. At launch, METAHICT detects the CPUs
and physical or cgroup-limited memory available to the local executor. Nextflow
uses those values as `resourceLimits`: any larger YAML or command-line request
is capped instead of stopping the workflow. For example, a 64 GB request on a
host with 48 GB becomes an effective 48 GB task allocation and produces a
warning.

The configured value remains the desired allocation; `task.cpus`,
`task.memory`, tool arguments, and the trace use the effective capped value. A
cap permits the task to start, but it cannot guarantee that a
memory-intensive method will succeed with less RAM than recommended.

## What threads and memory actually control

`resources.<stage>.threads` has two roles:

1. Nextflow requests that many CPUs and caps the request to local capacity; and
2. METAHICT passes the effective thread count to parallel tools in that stage.

For example, the binning value reaches binners and CheckM2, the annotation
value reaches GTDB-Tk, and the MGE value reaches geNomad and ccfind. A program
that is inherently serial cannot use extra CPUs; contact normalization
requests one.

`resources.<stage>.memory` is the requested Nextflow task allocation. Where a
tool exposes a meaningful total-memory flag, METAHICT also derives an
internal limit from 80% of the effective allocation after local capping:

- BBTools in preprocessing;
- MEGAHIT or metaSPAdes in assembly;
- BBMap in coverage; and
- SPAdes in reassembly.

The remaining 20% covers the process runtime and companion commands. Tools
such as BWA, CheckM2, YaHS, GTDB-Tk, geNomad, and ccfind do not share one
reliable total-memory flag. For them, the YAML memory is an executor request.

CPU and memory allocations for each task are recorded in
`results/nextflow_reports/trace.txt`.

## Resource precedence

The order is:

```text
command line -t/-m  >  YAML resources.<stage>  >  built-in fallback
```

For one selected stage, a command-line override is convenient:

```bash
./metahict run \
  --entry-module assembly \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --preprocessing-dir results/sample_01/1_preprocessing \
  --outdir results/assembly \
  -t 24 \
  -m 96G
```

For a complete workflow, `-t` and `-m` override every selected stage. Thus
`-m 96G` gives 96 GB to both light and heavy tasks, replacing all ten YAML
memory values. Per-stage YAML editing is normally better for a complete run.

The YAML omits wall-time limits because they depend on scheduler and site
policy. Scheduler users can define them in an institution-specific Nextflow
configuration.

## Method-specific settings

Only the subsection for the selected method is active.

### Assembly

```yaml
modules:
  assembly:
    assembler: megahit
    megahit:
      min_kmer_length: 21
    metaspades:
      kmer_lengths: "21,33,55"
    metaflye:
      extra_args: ""
```

With `assembler: megahit`, the `metaspades` map is retained as a documented
short-read alternative but is not passed to MEGAHIT. A samplesheet row with
`long_read_type` automatically selects metaFlye and uses `metaflye.extra_args`.

### Contact normalization

```yaml
modules:
  contact:
    method: normcc
    normcc:
      epsilon: 1
    hiczin:
      coverage_file: null
```

Valid methods are `raw`, `normcc`, `hiczin`, `bin3c`, and `metator`.
HiCzin and MetaTOR require a coverage table whose contig identifiers match
the selected assembly. The raw-contact and denoising maps apply before the
method-specific step.

### Binning

`modules.binning.metacc`, `.bin3c`, and `.imputecc` hold method-specific
settings. `refinement` controls quality evaluation and consolidation across
their candidate bins. `output` controls reports and optional FASTA products.
Changing a setting for one binner does not change the others.

### Reassembly read selection

The reassembly EM model treats the lower insert-size component as N and the
upper component as C. Its initialization and stopping criteria are configured
under `modules.reassembly.read_selection`:

```yaml
modules:
  reassembly:
    read_selection:
      em_initial_n_fraction: 0.8
      em_convergence_tolerance: 0.01
      em_max_iterations: 100
```

With these defaults, the lower 80% of mapped insert sizes initializes N and
the upper 20% initializes C. Fitting stops when the absolute change in
log-likelihood is below `0.01`, or after `100` iterations. The effective
values are recorded in `7_reassembly/run_parameters.yaml` and
`7_reassembly/read_selection_summary.json`.

### MGE

All MGE settings are under one `modules.mge` map and one `resources.mge` row.
The top-level `genomad`, `ccfind`, and `pairs` maps control MGE detection,
circular-contig discovery, and MGE–host pairing. Nested `alignment` and
`contact` maps are used only when METAHICT must generate its own Hi-C evidence.
By default, their BAM and contact matrices remain in the Nextflow work cache
and are not published as final results. Set `modules.mge.keep_intermediates:
true` to copy them under `10_MGE/intermediates/`.

## Preprocessing library switches

Shotgun and Hi-C reads share trimming and QC settings, while four booleans
control library-specific behavior:

```yaml
modules:
  preprocessing:
    libraries:
      shotgun:
        enabled: true
        deduplicate: false
      hic:
        enabled: true
        deduplicate: true
```

Defaults process both libraries, do not deduplicate short-read shotgun reads,
and do deduplicate Hi-C reads. A selected preprocessing run can disable one
library. A complete run keeps both library switches enabled; long-read rows
automatically bypass shotgun preprocessing.

## What belongs in the samplesheet

The minimal header is:

```csv
sample,sg1,sg2,hic1,hic2,enzyme,long_read_type
```

For paired short reads, `sg1` and `sg2` are required and `long_read_type` is
empty. For long reads, `sg1` is the single FASTA/FASTQ file, `sg2` is empty,
and `long_read_type` is one of `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`,
`nano-raw`, `nano-corr`, or `nano-hq`.

The enzyme is a sample/library fact and belongs in the CSV. It is
required by contact generation, binning, scaffolding, and analysis of
MGE–host pairs.
Algorithm parameters belong in YAML; the CSV contains no `*_extra_args`
columns.

## Module parameter reference

Each module page lists all keys, defaults, applicability, and output:

```bash
./metahict run --entry-module preprocessing --help
./metahict run --entry-module contact --help
./metahict run --entry-module binning --help
```

The same text is maintained under [Module reference](modules/README.md).
