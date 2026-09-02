# Module reference

The complete workflow runs nine core modules in dependency order. Scaffolding
is the tenth available module and is run separately on a selected MAG.
Selected-module runs use the same Nextflow resource management, exact
environments, reports, logs, and resume behavior as the complete workflow.

Definitions and workflow context are provided in [Concepts for new
users](../concepts.md).

## Stages at a glance

| Entry name | Purpose | Primary input | Guide |
| --- | --- | --- | --- |
| `preprocessing` | Clean and quality-check short reads | Raw short-read shotgun and paired Hi-C FASTQs | [Preprocessing](preprocessing.md) |
| `assembly` | Assemble contigs | Cleaned paired short reads or one long-read file | [Assembly](assembly.md) |
| `alignment` | Align Hi-C reads to contigs | Assembly and cleaned Hi-C reads | [Alignment](alignment.md) |
| `coverage` | Measure shotgun depth | Assembly and paired short reads or one long-read file | [Coverage](coverage.md) |
| `contact` | Construct and normalize Hi-C contacts | Assembly, Hi-C BAM, and enzyme | [Contact](contact.md) |
| `binning` | Recover and consolidate MAGs | Assembly, Hi-C BAM, enzyme, and quality databases | [Binning](binning.md) |
| `reassembly` | Recruit reads and improve short-read MAG assemblies | Bins and stages 1–3; paired short reads only | [Reassembly](reassembly.md) |
| `scaffolding` | Order/orient contigs in one MAG | One bin FASTA and matching Hi-C evidence | [Scaffolding](scaffolding.md) |
| `annotation` | Assign GTDB taxonomy | Directory of genome FASTAs | [Annotation](annotation.md) |
| `mge` | Detect MGEs, find circular contigs, and report candidate MGE–host pairs | Metagenome FASTA, host MAGs, and Hi-C evidence | [MGE](mge.md) |

## Get detailed terminal help

```bash
./metahict run --entry-module preprocessing --help
./metahict run --entry-module binning --help
```

For each entry, terminal help is generated from the corresponding maintained
manual and includes:

- what the stage does;
- required command-line inputs;
- every `modules.<stage>` YAML key and default;
- resource behavior;
- published outputs; and
- a complete selected-stage command.

## Common selected-module form

```bash
./metahict run \
  --entry-module NAME \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results/NAME \
  REQUIRED_UPSTREAM_OPTIONS
```

Replace `NAME` and `REQUIRED_UPSTREAM_OPTIONS` with the values listed in the
stage guide.

The samplesheet provides sample identity and, for contact-aware stages, the
enzyme. Each stage consumes only the relevant columns. Annotation, for
example, uses the sample name together with `--mag-dir`.

## Where settings belong

| Kind of information | Location |
| --- | --- |
| Sample name, read paths, Hi-C enzyme | `samplesheet.csv` |
| Normal per-stage threads and memory | `resources.<stage>` in YAML |
| Scientific and method-specific parameters | `modules.<stage>` in YAML |
| Upstream result selected for this invocation | `./metahict run` options |
| Database outside the default installation | `./metahict run --*-db` option |
| Temporary run-wide resource override | `-t/--threads`, `-m/--memory` |

Command-line resource values override YAML values for every selected task in
that invocation. For normal complete runs, configure each stage independently
in YAML.

## Interpreting success

A `[PASS]` workflow message means all requested processes exited successfully
and their declared outputs were published. Scientific acceptance also depends
on the stage-specific QC described in each guide. Retain the run record with
the reviewed results.

The public interface uses the descriptive entry names shown above with
`./metahict run`.
