# Testing METAHICT

METAHICT provides three levels of validation. Most users need the first two;
the complete example is recommended before analyzing important datasets.

| Test | Command | What it checks |
| --- | --- | --- |
| Installed runtime | `./metahict verify` | Exact Conda packages and pinned artifacts |
| Workflow stub | `./metahict test workflow` | Core workflow, standalone scaffolding entry, and output structure |
| Bundled example | `./metahict test example --outdir results` | Real core workflow followed by optional scaffolding |

`./metahict test source` is a developer test for source code, interfaces,
documentation, and policies. It is not required after an ordinary
installation.

## 1. Test the installation

Run immediately after `./metahict install`:

```bash
./metahict verify
./metahict test workflow
```

The workflow test generates temporary synthetic reads and dummy database
paths. Nextflow first runs the default core workflow in stub mode and then
runs the standalone scaffolding entry on the generated stub MAGs. It validates
the published files without executing the bioinformatics programs, so the
reference databases are not needed.

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

Run the core workflow, scaffold the recovered MAGs, and validate the outputs:

```bash
./metahict test example --outdir results
```

This command runs the actual scientific programs and databases. It first runs
the normal complete workflow, which excludes scaffolding, and then invokes the
standalone scaffolding entry for every MAG recovered from the example. It can
require substantial time, memory, and temporary storage even though the FASTQ
files are small.

The main results are written under:

```text
results/example_dataset/
```

The directory normally contains all ten numbered module directories because
the example test adds optional scaffolding after the core run. The number of
MAGs is not fixed. An unsuitable MAG is recorded as skipped instead of failing
the test; if no MAG is available, the command reports that scaffolding was not
run. Logs and provenance are under `results/nextflow_reports/`. A successful
test ends with:

```text
[PASS] Bundled example test completed
```

The validation checks the core output contract and the status of each attempted
scaffolding run. It accepts both completed and biologically ineligible
scaffolding outcomes. This confirms functional execution; it is not a
scientific accuracy benchmark.

## If a test fails

For a workflow or example failure, inspect:

```bash
tail -n 150 results/nextflow_reports/run.log
cat results/nextflow_reports/failure_summary.txt
```

The stub test uses a temporary directory, so its failing work path is printed
directly in the terminal. For the example run, follow the instructions in
[Troubleshooting](troubleshooting.md). After correcting the problem, keep
`results/nextflow_work/` and resume with:

```bash
./metahict test example --outdir results --resume
```

## Developer validation

After changing source code, workflow interfaces, configuration contracts,
tests, or manuals, run:

```bash
./metahict test source
./metahict test workflow
```

Changes that may affect scientific results should also run the bundled example
and compare its results with an accepted baseline using `./metahict compare`.
