# Running the workflow

METAHICT uses Nextflow to connect stages, reserve resources, isolate
dependencies, publish outputs, cache successful work, and generate reports.
`./metahict run` is the command-line interface for complete and selected-stage
analyses.

## Complete workflow

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results
```

The default entry is `all`. The launcher validates the configuration,
required database paths, and exact environments before Nextflow starts.

## Output and work directories

`--outdir` is the publication root. With sample `sample_01`, the structure is:

```text
results/
├── sample_01/
│   ├── 1_preprocessing/
│   ├── 2_assembly/
│   ├── 3_alignment/
│   ├── 4_coverage/
│   ├── 5_contact/
│   ├── 6_binning/
│   ├── 7_reassembly/
│   ├── 8_scaffolding/
│   ├── 9_annotation/
│   └── 10_MGE/
├── nextflow_reports/
└── nextflow_work/
```

Published scientific results and work files have different roles:

- keep the numbered sample directories as analysis results;
- keep `nextflow_reports/runs/<RUN_ID>/` as provenance;
- keep `nextflow_work/` while the run may need `--resume`;
- remove work files only after results and provenance have been validated and
  no resume is needed.

For long-read shotgun samples, stage 1 contains only `hic/` and stage 7 is
absent. The complete workflow uses the original long reads for stages 2 and 4,
then supplies stage 6 MAGs directly to stages 8–10.


## Run one module

Use a descriptive entry name and provide its upstream results. The table lists
the input contract for each entry.

| Entry | Required input beyond samplesheet/config | Typical source |
| --- | --- | --- |
| `preprocessing` | Raw reads already listed in samplesheet | Sequencing files |
| `assembly` | `--preprocessing-dir` for short reads; no upstream directory for long reads | Short reads use `1_preprocessing/sg`; long reads use samplesheet `sg1` |
| `alignment` | `--assembly-dir`, `--preprocessing-dir` | `2_assembly`, `1_preprocessing`; alignment uses `hic/` |
| `coverage` | `--assembly-dir`; also `--preprocessing-dir` for short reads | Coverage maps cleaned short reads or samplesheet long reads |
| `contact` | `--assembly-dir`, `--alignment-dir`; enzyme in samplesheet | `2_assembly`, `3_alignment` |
| `binning` | `--assembly-dir`, `--alignment-dir`; enzyme and databases | `2_assembly`, `3_alignment` |
| `reassembly` | `--binning-dir`, `--assembly-dir`, `--alignment-dir`, `--preprocessing-dir` | Paired short-read samples only; stages 1, 2, 3, and 6 are used |
| `scaffolding` | `--scaffolding-bin`, `--preprocessing-dir`; enzyme | One MAG FASTA and the `hic/` child of preprocessing |
| `annotation` | `--mag-dir` | Directory directly containing MAG FASTAs |
| `mge` | `--fasta`, `--host-dir`, and either cleaned Hi-C reads or reusable MGE contact/alignment results | Matching sequences, host genomes, and Hi-C evidence for MGE detection, circular-contig discovery, and MGE–host pairs |

Example alignment run:

```bash
./metahict run \
  --entry-module alignment \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --preprocessing-dir results/sample_01/1_preprocessing \
  --assembly-dir results/sample_01/2_assembly \
  --outdir results/alignment
```

Get the full command, inputs, parameters, resource behavior, and outputs for
any stage:

```bash
./metahict run --entry-module reassembly --help
```


## Special input rules

### Scaffolding accepts one bin

The selected scaffolding entry operates on the FASTA supplied with
`--scaffolding-bin`. It aligns cleaned Hi-C reads to that bin unless
`--scaffolding-bam` supplies a BAM aligned to exactly the same reference
sequences and lengths. The complete workflow automatically creates one
scaffolding task for every reassembled bin.

### Annotation accepts a MAG directory

`--mag-dir` must directly contain `.fa`, `.fasta`, or `.fna` MAG files.
Annotation needs a samplesheet for sample identity, but it does not consume
read data or infer an extra `reassembled_bins` directory.

For multiple samples, `{sample}` is replaced with each sample identifier:

```bash
./metahict run \
  --entry-module annotation \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --mag-dir results/reassembly/{sample}/7_reassembly/reassembled_bins \
  --outdir results/annotation
```

### MGE accepts generic metagenome FASTAs

`--fasta` can be a short-read, long-read, hybrid, or reassembled metagenome.
`--host-dir` supplies the matching host MAGs, and Hi-C evidence must use the
same contig identifiers. Every run detects MGEs, finds circular contigs, and
reports candidate MGE–host pairs.

## Resource overrides

Normal per-stage resources belong in `metahict_configuration.yaml`. A one-run
`-t` or `-m` value has highest precedence and replaces every selected stage's
configured value:

```bash
./metahict run \
  --entry-module binning \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --assembly-dir results/sample_01/2_assembly \
  --alignment-dir results/sample_01/3_alignment \
  --outdir results/binning \
  -t 24 \
  -m 112G
```

For a complete run, edit the relevant YAML row so that each stage retains an
appropriate allocation.

On the supported local executor, `./metahict` detects the CPUs and memory
available to the process and supplies them as Nextflow `resourceLimits`. If a
configured or command-line request is larger, Nextflow warns, caps the task,
and passes the effective `task.cpus` and `task.memory` values to the module.

## Preview without running

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results \
  --show-command
```

This prints the generated Nextflow command after path resolution and exits.
It checks command construction; compatibility between independently generated
scientific files still requires review.

## Resume

After correcting an error, repeat the same inputs, configuration, output root,
and work directory with `--resume`:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results \
  --resume
```

Nextflow reuses compatible successful tasks from `nextflow_work/` and reruns
failed or invalidated work. Changing an input, parameter, script, or command
may invalidate downstream cache entries.

## Run outside the source directory

METAHICT can be started from any directory. Use the launcher's absolute path
and absolute analysis paths:

```bash
/path/to/METAHICT/metahict run \
  --samplesheet /path/to/analysis/samplesheet.csv \
  --config /path/to/analysis/metahict_configuration.yaml \
  --outdir /path/to/analysis/results
```

The launcher locates its own workflow and environments. Relative command-line
paths are resolved from the current directory, so absolute paths are safest
for saved server commands.

## Common run options

| Option | Meaning |
| --- | --- |
| `--samplesheet` | Sample identities, raw-read paths, and enzymes |
| `--config` | All-module YAML resources and scientific settings |
| `--entry-module` | `all` or one descriptive stage name |
| `--outdir` | Published result root |
| `--report-dir` | Override the immutable report root |
| `--work-dir` | Override the Nextflow cache/work root |
| `-t`, `--threads` | Run-wide thread override |
| `-m`, `--memory` | Run-wide memory override, such as `112G` |
| `--resume` | Reuse compatible successful tasks |
| `--verbose-preflight` | Print every environment verification line |
| `--show-command` | Print the generated Nextflow command and exit |

Run `./metahict run --help` for database, scaffolding, annotation, and MGE
options.

## Direct Nextflow use

Workflow developers can inspect `nextflow/main_dsl2.nf`,
`nextflow/modules/local/metahict_modules.nf`, and the command emitted by
`--show-command`. Direct Nextflow invocation bypasses parts of CLI validation
and run-record management.
