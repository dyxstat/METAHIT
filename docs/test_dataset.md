# Testing METAHICT

METAHICT provides three levels of validation. Most users need the first two;
the complete example is recommended before analyzing important datasets.

| Test | Command | What it checks |
| --- | --- | --- |
| Installed runtime | `./metahict verify` | Exact Conda packages and pinned artifacts |
| Workflow stub | `./metahict test workflow` | Nextflow compilation, all stage connections, and output structure |
| Bundled example | `./metahict run ... --check-outputs` | Real execution of every scientific stage and required database |

`./metahict test source` is a developer test for source code, interfaces,
documentation, and policies. It is not required after an ordinary
installation.

## 1. Test the installation

Run immediately after `./metahict install`:

```bash
./metahict verify
./metahict test workflow
```

The workflow test generates temporary synthetic paired reads and dummy
database paths. Nextflow runs every process in stub mode and validates the
published files. It does not execute the bioinformatics programs on real data,
so the reference databases are not needed.

The test passes when both Nextflow stub profiles and the expected-output check
end with `[PASS]`.

## 2. Install and check the databases

```bash
./metahict database all
./metahict doctor --runtime --databases
```

Resolve every failed environment or database check before running the bundled
example.

## 3. Run the bundled example

The repository includes synchronized paired shotgun and Hi-C reads under
`example_dataset/`. They are a compact functional test, not an unbiased
biological benchmark.

Run every stage and validate the documented outputs:

```bash
./metahict run \
  --samplesheet nextflow/assets/example_dataset_samplesheet.csv \
  --config nextflow/assets/example_dataset_configuration.yaml \
  --outdir results \
  --check-outputs
```

This command runs the actual scientific programs and databases. It can require
substantial time, memory, and temporary storage even though the FASTQ files are
small.

The main results are written under:

```text
results/example_dataset/
```

The directory should contain all ten numbered stages. Logs and provenance are
under `results/nextflow_reports/`. A successful run ends with:

```text
[PASS] METAHICT workflow completed
```

`--check-outputs` then verifies that the required files exist and are
non-empty. This confirms functional execution; it is not a scientific accuracy
benchmark.

## If a test fails

For a workflow or example failure, inspect:

```bash
tail -n 150 results/nextflow_reports/run.log
cat results/nextflow_reports/failure_summary.txt
```

The stub test uses a temporary directory, so its failing work path is printed
directly in the terminal. For the example run, follow the instructions in
[Troubleshooting](troubleshooting.md). After correcting the problem, repeat the
example command with `--resume` and keep `results/nextflow_work/`.

## Developer validation

After changing source code, workflow interfaces, configuration contracts,
tests, or manuals, run:

```bash
./metahict test source
./metahict test workflow
```

Changes that may affect scientific results should also run the bundled example
and compare its results with an accepted baseline using `./metahict compare`.
