# Troubleshooting

Start with the first failed process rather than the last repeated Python
exception. METAHICT preserves the information needed to diagnose it.

## First five checks

For the output root `results`, run:

```bash
tail -n 80 results/nextflow_reports/run.log
cat results/nextflow_reports/failure_summary.txt
grep -E 'FAILED|ABORTED' results/nextflow_reports/trace.txt
head -n 20 results/nextflow_reports/process_logs/index.tsv
./metahict doctor --runtime --databases
```

The failure summary identifies the archived log directory when Nextflow
recorded the failed task. Inspect `command.err`, `command.out`, and
`command.sh` there. See [Logging and provenance](logging.md) for the complete
layout.

## A resource request is larger than the machine

Example:

```text
[WARN] Requested resources exceed local capacity: up to 16 threads/64 GB requested; 8 threads/48 GB detected. Nextflow will cap each task to the detected limits and pass the effective values to its tools.
```

This warning is informational. METAHICT passes the detected CPU and memory
ceilings to Nextflow's local-executor `resourceLimits`. Nextflow reduces an
oversized task request, starts the task, and exposes the effective values as
`task.cpus` and `task.memory`; METAHICT then passes those values to parallel and
memory-aware tools.

The task can still fail if the underlying method genuinely needs more memory.
GTDB-Tk full-tree classification, large assemblies, and large contact matrices
are examples where a larger machine may still be necessary. Check the effective
allocation in `results/nextflow_reports/trace.txt` if a capped task fails.

## The process was killed or exited with code 137

This commonly indicates memory exhaustion. Compare the task's `memory` in
`trace.txt` with host availability and inspect `command.err`. Increase the
configured memory, reduce tool concurrency where supported, or use a larger
host. Some programs allocate memory per thread, so increasing threads can also
increase memory use.

## A FASTQ or upstream directory cannot be found

Use `./metahict samplesheet` so FASTQ paths are written as absolute paths.
Confirm files and gzip integrity:

```bash
ls -lh /data/sample_01/shotgun_R1.fastq.gz
gzip -t /data/sample_01/shotgun_R1.fastq.gz
```

For a selected module, pass the directory that directly contains the expected
published files. New outputs are flat; for example, assembly input is
`results/sample_01/2_assembly`, not an extra
`2_assembly/assembly` level. The module-specific help shows every required
directory:

```bash
./metahict run --entry-module reassembly --help
```

## The restriction enzyme is missing

Contact, binning, scaffolding, and analysis of MGE–host pairs require the `enzyme` column
in the samplesheet. Regenerate or correct the row using the library
preparation record. There is no scientifically safe default enzyme.

## An environment lock does not match

Run:

```bash
./metahict verify
```

The detailed output names the mismatched project-local environment. Preserve
anything needed from that directory, move or remove only the named
environment, then rerun `./metahict install`. `metahict_env` is the main
scientific runtime and should be repaired through the installer.

## BBTools cannot find `adapters.fa`

Run `./metahict verify`. The installer verifies the Conda-installed BBTools
adapter reference and the workflow resolves it from the installed package.
If verification fails, repair the `metahict_env` installation. A custom
adapter reference belongs at
`modules.preprocessing.trimming.adapter_reference` in the YAML file.

## CheckM2 says DIAMOND is missing

Verify the executable in the CheckM2 environment:

```bash
conda run -p conda_envs/checkm2 diamond version
./metahict verify
```

If DIAMOND is absent, repair the locked `checkm2` environment. METAHICT adds
that environment to the task `PATH`; installing an unrelated system DIAMOND
does not repair an incomplete lock installation.

## ccfind says `ssearch36` is missing

Run `./metahict verify`. `ssearch36` is supplied by the locked `ccfind_env`,
which METAHICT adds to the MGE task `PATH`. Repair that environment if the
check fails; the installer restores the expected executable and version.


## Resume after fixing the problem

Repeat the same command with the same output and work directories and add
`--resume`:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results \
  --resume
```

Preserve `results/nextflow_work/`; it contains the resume cache.
Nextflow reruns failed or invalidated tasks and reuses compatible successful
tasks.

## Information to include in a support request

Include:

- `./metahict --version`;
- operating system and available CPUs and RAM;
- `provenance.json`, `run.log`, and `failure_summary.txt`;
- the failed row from `trace.txt`; and
- the failed task's archived `command.sh`, `command.err`, and `command.out`.

Redact confidential reads, credentials, and identifiable host paths from
support files.
