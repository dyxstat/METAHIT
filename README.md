# METAHICT

METAHICT is a Nextflow DSL2 workflow for genome-resolved analysis of short- or
long-read shotgun metagenomes with paired-end metagenomic Hi-C data. It performs read preprocessing,
assembly, Hi-C alignment, coverage and contact analysis, binning, bin
reassembly, scaffolding, taxonomic annotation, and detection of mobile genetic
elements (MGEs) and candidate MGE–host pairs.

![METAHICT workflow overview](images/METAHICT_Overview.png)

## Quick start

### Requirements

METAHICT 1.1.0 supports and tests on 64-bit x86 Linux systems. The host must
provide Conda, `curl`, `tar`, and Git. Python, Nextflow, and the scientific
programs are installed into project-local environments by METAHICT.

Default allocations reach 16 threads and 64 GB RAM. METAHICT can cap requests
to the available local resources, but large datasets may still need a
high-memory server. Reference databases, results, and Nextflow work files also
require substantial storage.

### 1. Download METAHICT

```bash
git clone https://github.com/dyxstat/METAHICT.git
cd METAHICT
chmod +x metahict nextflow/bin/nextflow
```

### 2. Install and test the software

```bash
./metahict doctor
./metahict install
./metahict test workflow
```

`doctor` checks the operating system, architecture, required commands, and
distributed lock files before installation. `install` creates and verifies the
locked software environments. `test workflow` runs the core workflow and the
standalone scaffolding entry with stub tasks, without downloading databases or
analyzing real reads. Typical installation time is 1–5 minutes, excluding
database installation.


### 3. Install and validate the reference databases

```bash
./metahict database all
./metahict doctor --runtime --databases
```

This installs the CheckM, CheckM2, GTDB-Tk, and geNomad databases in the
default `databases/` directory. Shared or existing databases are supported;
see [Installation and databases](docs/installation.md).

### 4. Create the configuration and samplesheet

Create the complete configuration template:

```bash
./metahict config
```

For a first analysis, keep the scientific defaults and review the `resources`
section of `metahict_configuration.yaml` for the available server.

Create a samplesheet using the real FASTQ paths and the restriction enzyme or
enzyme combination used to prepare the Hi-C library:

```bash
./metahict samplesheet \
  --sample sample_01 \
  --sg-r1 /data/sample_01/shotgun_R1.fastq.gz \
  --sg-r2 /data/sample_01/shotgun_R2.fastq.gz \
  --hic-r1 /data/sample_01/hic_R1.fastq.gz \
  --hic-r2 /data/sample_01/hic_R2.fastq.gz \
  --enzyme Sau3AI,MluCI \
  --output samplesheet.csv
```

Replace the paths, sample name, and enzymes with the study metadata.

For a single-file long-read shotgun library, omit `--sg-r2` and add its
long read type, for example `--long-read-type nano-hq`. Accepted values are
`pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, and
`nano-hq`.

### 5. Run the complete workflow

Start the analysis:

```bash
./metahict run \
  --samplesheet samplesheet.csv \
  --config metahict_configuration.yaml \
  --outdir results
```

METAHICT writes each sample to `results/<sample>/` and workflow reports to
`results/nextflow_reports/`. Scaffolding is not run automatically; it remains
available as a selected-module analysis for MAGs chosen by the user.

After correcting a failed task, repeat the run command with `--resume` and
keep `results/nextflow_work/` so Nextflow can reuse successful tasks.

### 6. Run individual modules (optional)

METAHICT can run one stage independently when its required upstream results
already exist. Available entry names are `preprocessing`, `assembly`,
`alignment`, `coverage`, `contact`, `binning`, `reassembly`, `scaffolding`,
`annotation`, and `mge`.

Before running a stage, display its required inputs, configuration parameters,
outputs, and complete example command:

```bash
./metahict run --entry-module binning --help
```

Replace `binning` with the required entry name. The [module
reference](docs/modules/README.md) provides the same information as browsable
guides for all ten stages.

## Test with the bundled example dataset

After installing the databases, run the core workflow and then exercise
scaffolding on every MAG recovered from the included paired reads:

```bash
./metahict test example --outdir results
```

This is a real scientific-program test and can take substantially longer than
the workflow stub. The bundled example test typically takes approximately
25–30 minutes, depending on the system. Details are in [Testing
METAHICT](docs/test_dataset.md).

## Main results

| Result | Location below `results/sample_01/` |
| --- | --- |
| Cleaned shotgun and Hi-C reads | `1_preprocessing/sg/`; `1_preprocessing/hic/` |
| Metagenome assembly | `2_assembly/final_assembly.fasta` |
| Filtered Hi-C alignment | `3_alignment/sorted_map.bam` |
| Shotgun depth | `4_coverage/coverage.txt` |
| Normalized contact matrix | `5_contact/denoised_contact_matrix_normcc.npz` |
| Consolidated MAGs | `6_binning/metahict/final_bins/` |
| Reassembled MAGs (paired short reads only) | `7_reassembly/reassembled_bins/` |
| Scaffolded MAGs (optional standalone module) | `8_scaffolding/<BIN_ID>/scaffolded_bin.fa` |
| GTDB-Tk taxonomy | `9_annotation/classify/gtdbtk.*.summary.tsv` |
| MGE calls, circular contigs, and candidate MGE–host pairs | `10_MGE/` |

## Run from another directory

The launcher can be called by its absolute path. Absolute analysis paths make
the saved command unambiguous:

```bash
/path/to/METAHICT/metahict run \
  --samplesheet /path/to/analysis/samplesheet.csv \
  --config /path/to/analysis/metahict_configuration.yaml \
  --outdir /path/to/analysis/results
```

## Documentation and help

| Goal | Documentation |
| --- | --- |
| First complete analysis | [Command-by-command tutorial](docs/quickstart.md) |
| Understand the biological stages | [Concepts](docs/concepts.md) |
| Install or reuse environments and databases | [Installation and databases](docs/installation.md) |
| Change resources or algorithms | [Configuration reference](docs/configuration.md) |
| Run one module | [Module reference](docs/modules/README.md) |
| Understand results, resume, or run outside the checkout | [Workflow execution](docs/nextflow.md) |
| Inspect logs and provenance | [Logging](docs/logging.md) |
| Diagnose a failure | [Troubleshooting](docs/troubleshooting.md) |

Command-line help is generated from the current interface:

```bash
./metahict --help
./metahict run --help
./metahict run --entry-module binning --help
```

The complete documentation index is [docs/README.md](docs/README.md).

## Third-party software

Versions, licenses, sources, and citations for external programs are listed in
[Third-party software](docs/third_party.md).

## License

METAHICT is distributed under the GNU General Public License; see
[LICENSE](LICENSE). A workflow-layer Bioconda recipe is staged under
`packaging/bioconda/` but is not a published Bioconda package until accepted by
the Bioconda project.
