# METAHICT documentation

The shortest path from download to a complete analysis is:

1. [Command-by-command tutorial](quickstart.md)
2. [Configuration reference](configuration.md)
3. [Troubleshooting](troubleshooting.md), if a task fails

The repository [README](../README.md) contains a compact installation and run
example. The pages below provide additional explanation when it is needed.

## Getting started

| Question | Guide |
| --- | --- |
| What information do shotgun and Hi-C reads provide? | [Concepts](concepts.md) |
| How do I install METAHICT or reuse an existing installation? | [Installation and databases](installation.md) |
| What commands do I run for my first sample? | [Command-by-command tutorial](quickstart.md) |
| How do I confirm that the installation works? | [Testing METAHICT](test_dataset.md) |

## Running an analysis

| Question | Guide |
| --- | --- |
| How do I set threads, memory, and scientific methods? | [Configuration](configuration.md) |
| How do I run the complete workflow or one module? | [Workflow execution](nextflow.md) |
| What are the inputs and outputs of each module? | [Module reference](modules/README.md) |
| Where are logs, reports, and provenance records? | [Logging](logging.md) |
| How do I diagnose and resume a failed run? | [Troubleshooting](troubleshooting.md) |

## Technical reference

- [Architecture](architecture.md): responsibilities of the launcher,
  Nextflow, Python stages, and external programs.
- [Third-party software](third_party.md): versions, licenses, sources,
  installation routes, and citations.

The current command-line interface is also a reference:

```bash
./metahict --help
./metahict run --help
./metahict run --entry-module MODULE --help
```

Replace `MODULE` with `preprocessing`, `assembly`, `alignment`, `coverage`,
`contact`, `binning`, `reassembly`, `scaffolding`, `annotation`, or `mge`.
