# METAHICT

**METAHICT** enables comprehensive genome-resolved microbiome analysis with metagenomic Hi-C.

METAHICT v1.1.0 is distributed as a Nextflow DSL2 workflow with module-level command-line wrappers. The workflow supports raw-read processing, assembly, Hi-C alignment, coverage estimation, contact generation and normalization, Hi-C-informed binning and consolidation, per-bin reassembly, Hi-C-guided scaffolding, GTDB-Tk annotation, and MGE-MAG candidate proximity-association analysis.

![METAHICT overview](images/METAHICT_Overview.png)

## Quick start

Clone the repository:

```bash
git clone https://github.com/dyxstat/METAHICT.git
cd METAHICT
```

Install the local software environments:

```bash
bash installation/run_setup_in_venv.sh
```

Download the reference databases into the default project layout:

```bash
bash installation/db/checkm_db.sh databases/checkm_db
bash installation/db/checkm2_db.sh databases/checkm2_db
bash installation/db/gtdbtk_db.sh databases/gtdbtk_db/release220
bash installation/db/genomad_db.sh databases/genomad_db
```

Run the bundled Nextflow smoke test:

```bash
bash nextflow/ci/run_smoke_ci.sh
```

After downloading the example FASTQ files from Zenodo into `example_data/`, run the workflow with the example sample sheet:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/example_data" \
  --report_dir "$PWD/results/example_data/nextflow_reports" \
  -work-dir "$PWD/results/example_data/nextflow_work" \
  -ansi-log false
```

## Installation

### Local Conda installation

The local installation is intended for systems where containers are not available:

```bash
bash installation/run_setup_in_venv.sh
```

This creates a project-local bootstrap environment and installs the METAHICT environments under `conda_envs/`. 

The local install creates these environments:

| Environment | Main use |
| --- | --- |
| `metahict_env` | Core METAHICT wrappers and most module tools |
| `checkm2` | CheckM2 quality evaluation |
| `gtdbtk-2.4.0` | GTDB-Tk annotation |
| `genomad` | geNomad MGE discovery |
| `ccfind_env` | ccfind circularity detection |

### Nextflow profiles

The release defines the following Nextflow profiles in `nextflow/nextflow.config`:

| Profile | Use case |
| --- | --- |
| `local` | Run with the project-local environments installed by `installation/run_setup_in_venv.sh` |
| `conda` | Verify and use the locked project-local Conda environment bundle |
| `docker` | Run processes with Docker |
| `apptainer` | Run processes with Apptainer, commonly used on HPC systems |
| `singularity` | Run processes with Singularity |
| `slurm` | Submit Nextflow tasks to a SLURM executor; requires site-specific cluster configuration |

For containerized runs, provide the METAHICT all-tools image with `--container_image_override`.

Example Apptainer command using the configured local SIF path:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile apptainer \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/example_data" \
  --report_dir "$PWD/results/example_data/nextflow_reports" \
  -work-dir "$PWD/results/example_data/nextflow_work" \
  --container_image_override "$PWD/nextflow/containers/images/metahict-all_local.sif" \
  -ansi-log false
```

## Databases

Large reference databases are not bundled in the Git repository. Download them separately and pass their paths through Nextflow parameters if they are not stored in the default project layout.

| Database | Download script | Used by |
| --- | --- | --- |
| CheckM | `installation/db/checkm_db.sh` | Hi-C binning through ImputeCC/CheckM marker profiling |
| CheckM2 | `installation/db/checkm2_db.sh` | Binning, reassembly, and scaffolding quality evaluation |
| GTDB-Tk release 220 | `installation/db/gtdbtk_db.sh` | Annotation |
| geNomad | `installation/db/genomad_db.sh` | MGE discovery and MGE taxonomy |

Database paths can be supplied to the workflow as:

```bash
--checkm_db /path/to/checkm_db
--checkm2_db /path/to/uniref100.KO.1.dmnd
--gtdbtk_db /path/to/gtdbtk/release220
--genomad_db /path/to/genomad_db
```

If the database scripts are run into the default `databases/` directory, the defaults in `nextflow/nextflow.config` can be used directly.

## Example dataset

The METAHICT example dataset is a 5% subset of the human-gut benchmark dataset. It contains paired-end shotgun reads and paired-end metagenomic Hi-C reads for installation checks, tutorial runs, and workflow smoke testing.

The example FASTQ files are archived at:

```text
https://doi.org/10.5281/zenodo.21695166
```

The repository includes the matching metadata and workflow input files:

| File | Purpose |
| --- | --- |
| `example_data/README.md` | Example-data description |
| `example_data/manifest.tsv` | Read counts for the example FASTQ files |
| `example_data/MD5SUMS.txt` | Checksums for downloaded example files |
| `nextflow/assets/example_data_samplesheet.csv` | Nextflow sample sheet for the example run |

The example sample sheet expects the FASTQ files to be placed under `example_data/` with these names:

```text
example_data/sg_R1.fastq.gz
example_data/sg_R2.fastq.gz
example_data/hic_R1.fastq.gz
example_data/hic_R2.fastq.gz
```

## Running METAHICT with Nextflow

The main workflow entry point is:

```text
nextflow/main_dsl2.nf
```

The wrapper script:

```text
nextflow/run_metahict_nextflow.sh
```

sets the project-specific Nextflow home, log location, Apptainer cache, and local PATH entries before launching the bundled Nextflow executable.

### Full workflow

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/example_data" \
  --report_dir "$PWD/results/example_data/nextflow_reports" \
  -work-dir "$PWD/results/example_data/nextflow_work" \
  -ansi-log false
```

### Per-module Nextflow execution

Use `--entry_module` to run one workflow module at a time. Valid values are `all` or `module1` through `module10`.

The commands below assume that previous module outputs are stored under:

```text
$PWD/results/example_data/example_data
```

Run module 1, preprocessing:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module1 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --out_root "$PWD/results/module1" \
  --report_dir "$PWD/results/module1/nextflow_reports" \
  -work-dir "$PWD/results/module1/nextflow_work" \
  -ansi-log false
```

Run module 2, assembly:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module2 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --sg_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/sg" \
  --out_root "$PWD/results/module2" \
  --report_dir "$PWD/results/module2/nextflow_reports" \
  -work-dir "$PWD/results/module2/nextflow_work" \
  -ansi-log false
```

Run module 3, Hi-C alignment:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module3 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --hic_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module3" \
  --report_dir "$PWD/results/module3/nextflow_reports" \
  -work-dir "$PWD/results/module3/nextflow_work" \
  -ansi-log false
```

Run module 4, shotgun coverage:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module4 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --sg_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/sg" \
  --out_root "$PWD/results/module4" \
  --report_dir "$PWD/results/module4/nextflow_reports" \
  -work-dir "$PWD/results/module4/nextflow_work" \
  -ansi-log false
```

Run module 5, contact generation and normalization:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module5 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --alignment_dir "$PWD/results/example_data/example_data/3_alignment" \
  --out_root "$PWD/results/module5" \
  --report_dir "$PWD/results/module5/nextflow_reports" \
  -work-dir "$PWD/results/module5/nextflow_work" \
  -ansi-log false
```

Run module 6, Hi-C binning and consolidation:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module6 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --alignment_dir "$PWD/results/example_data/example_data/3_alignment" \
  --out_root "$PWD/results/module6" \
  --report_dir "$PWD/results/module6/nextflow_reports" \
  -work-dir "$PWD/results/module6/nextflow_work" \
  -ansi-log false
```

Run module 7, reassembly:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module7 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --binning_dir "$PWD/results/example_data/example_data/6_binning" \
  --assembly_dir "$PWD/results/example_data/example_data/2_assembly" \
  --alignment_dir "$PWD/results/example_data/example_data/3_alignment" \
  --sg_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/sg" \
  --hic_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module7" \
  --report_dir "$PWD/results/module7/nextflow_reports" \
  -work-dir "$PWD/results/module7/nextflow_work" \
  -ansi-log false
```

Run module 8, scaffolding:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module8 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --reassembly_dir "$PWD/results/example_data/example_data/7_reassembly/reassembly" \
  --alignment_dir "$PWD/results/example_data/example_data/3_alignment" \
  --hic_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module8" \
  --report_dir "$PWD/results/module8/nextflow_reports" \
  -work-dir "$PWD/results/module8/nextflow_work" \
  -ansi-log false
```

Run module 9, annotation:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module9 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --reassembly_dir "$PWD/results/example_data/example_data/7_reassembly/reassembly" \
  --out_root "$PWD/results/module9" \
  --report_dir "$PWD/results/module9/nextflow_reports" \
  -work-dir "$PWD/results/module9/nextflow_work" \
  -ansi-log false
```

Run module 10, MGE analysis, including MGE-specific Hi-C alignment and contact generation:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module10 \
  --samplesheet nextflow/assets/example_data_samplesheet.csv \
  --reassembly_dir "$PWD/results/example_data/example_data/7_reassembly/reassembly" \
  --hic_preprocessing_dir "$PWD/results/example_data/example_data/1_preprocessing/hic" \
  --out_root "$PWD/results/module10" \
  --report_dir "$PWD/results/module10/nextflow_reports" \
  -work-dir "$PWD/results/module10/nextflow_work" \
  -ansi-log false
```

Run only the final module 10 MGE step when MGE-specific contact outputs already exist:

```bash
bash nextflow/run_metahict_nextflow.sh run nextflow/main_dsl2.nf \
  -profile local \
  --entry_module module10 \
  --samplesheet nextflow/assets/hg_mge_samplesheet.csv \
  --reassembly_dir "$PWD/hg/7_reassembly/results/reassembly_sg_hic" \
  --mge_contact_dir "$PWD/hg/10_MGE/results/mge_contact" \
  --out_root "$PWD/hg/10_MGE/rerun" \
  --report_dir "$PWD/hg/10_MGE/rerun/nextflow_reports" \
  -work-dir "$PWD/hg/10_MGE/rerun/nextflow_work" \
  -ansi-log false
```

Selected-module runs can reuse previous outputs through explicit directory parameters:

| Parameter | Meaning |
| --- | --- |
| `--sg_preprocessing_dir` | Existing shotgun preprocessing output |
| `--hic_preprocessing_dir` | Existing Hi-C preprocessing output |
| `--assembly_dir` | Existing assembly output |
| `--alignment_dir` | Existing Hi-C alignment output |
| `--coverage_dir` | Existing shotgun coverage output |
| `--binning_dir` | Existing binning output |
| `--reassembly_dir` | Existing reassembly output |
| `--mge_alignment_dir` | Existing MGE-specific Hi-C alignment output |
| `--mge_contact_dir` | Existing MGE-specific contact output |

### Sample sheet

Each Nextflow run uses a CSV sample sheet to define the sample name, input files, and Hi-C restriction enzyme information. The repository includes an example sample sheet for the Zenodo example dataset:

```text
nextflow/assets/example_data_samplesheet.csv
```

For raw short-read runs, the essential columns are:

| Column | Description |
| --- | --- |
| `sample` | Sample name |
| `sg1`, `sg2` | Shotgun paired-end FASTQ files |
| `hic1`, `hic2` | Hi-C paired-end FASTQ files |
| `enzyme` | Restriction enzyme list, for example `Sau3AI,MluCI` |

Users can copy this file and replace the paths with their own data. Per-module runs can also reuse existing output directories through command-line parameters such as `--assembly_dir`, `--alignment_dir`, `--binning_dir`, and `--reassembly_dir`.

### Output layout

For a sample named `example_data`, the full workflow publishes outputs under:

```text
<out_root>/example_data/
```

Major output directories include:

| Directory | Main contents |
| --- | --- |
| `1_preprocessing/sg/` and `1_preprocessing/hic/` | Cleaned shotgun and Hi-C reads |
| `2_assembly/` | Assembled contigs |
| `3_alignment/` | Hi-C BAM and 3D-ratio output |
| `4_coverage/` | Shotgun coverage BAM and coverage table |
| `5_contact/` | Raw and normalized Hi-C contact matrices |
| `6_binning/` | Individual binner outputs and METAHICT consolidated bins |
| `7_reassembly/` | Reassembled bins, name map, residual assembly, and combined contigs |
| `8_scaffolding/` | YaHS scaffolds, heatmap, and scaffolding metrics |
| `9_annotation/` | GTDB-Tk classification output |
| `10_MGE/` | MGE-specific alignment/contact outputs and MGE reports |

Nextflow provenance files are written to `--report_dir`, including `trace.txt`, `timeline.html`, `report.html`, and `dag.html`.

## Key outputs

### Consolidated MAGs

METAHICT binning outputs include:

```text
6_binning/binning/metahict/metahict_50_10_bins/
```

These FASTA files are the non-redundant consolidated MAG collection generated from the Hi-C-aware binning and refinement procedure.

### Reassembled bins

Reassembly outputs include:

```text
7_reassembly/reassembly/reassembled_bins/
7_reassembly/reassembly/reassembled_bin_name_map.tsv
7_reassembly/reassembly/combined_contigs.fa
```

The final reassembled bins use clean bin IDs, while `reassembled_bin_name_map.tsv` records whether the selected representative came from the original, strict, or permissive reassembly candidate.

### MGE-MAG candidate associations

The MGE module reports a single filtered candidate association table:

```text
10_MGE/mge/candidate_mge_mag_associations_zscore_filtered.tsv
```

The default filter retains candidate MGE-MAG proximity associations supported by at least two raw Hi-C read pairs and normalized-contact Z-score greater than 0.5. METAHICT reports these as candidate proximity associations, not confirmed infections or stable host ranges.

The MGE module also reports sequence topology for all combined contigs:

```text
10_MGE/mge/sequence_topology.tsv
```

This table has four columns:

| Column | Meaning |
| --- | --- |
| `contig_id` | Combined-contig identifier |
| `mge` | `1` if the contig is classified as viral or plasmid by geNomad, otherwise `0` |
| `ccfind` | `1` if ccfind reports circularity, otherwise `0` |
| `genomad` | `1` if geNomad reports circularity for an MGE contig, otherwise `0`; `NA` for non-MGE contigs |

The MGE module does not produce an unfiltered association table, a primary/best-host table, host-domain incompatibility columns, or a separate `circular_mge_candidates.tsv` file in v1.1.0.

## Manual module commands

The Nextflow workflow is the recommended user interface. The module-level command wrapper remains available for manual reruns, debugging, and selected analyses.

Print available commands:

```bash
python metahict.py --help
```

Print module defaults:

```bash
python metahict.py <module> --print-defaults
```

Per-module command templates:

```bash
python metahict.py preprocessing \
  -p "$PWD" -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o out/preprocessing -t 80

python metahict.py assembly \
  -p "$PWD" -1 out/preprocessing/final_sg_1.fastq.gz -2 out/preprocessing/final_sg_2.fastq.gz -o out/assembly -t 80

python metahict.py alignment \
  -p "$PWD" -r out/assembly/final_assembly.fasta \
  -1 out/preprocessing/final_hic_1.fastq.gz -2 out/preprocessing/final_hic_2.fastq.gz \
  -o out/alignment -t 80

python metahict.py coverage \
  -p "$PWD" -r out/assembly/final_assembly.fasta \
  -1 out/preprocessing/final_sg_1.fastq.gz -2 out/preprocessing/final_sg_2.fastq.gz \
  -o out/coverage -t 80

python metahict.py contact normcc \
  -p "$PWD" --bam out/alignment/sorted_map.bam --fasta out/assembly/final_assembly.fasta \
  --out out/contact --enzyme Sau3AI,MluCI

python metahict.py binning \
  out/assembly/final_assembly.fasta out/alignment/sorted_map.bam out/binning "$PWD" \
  -t 80 --checkm_db "$PWD/databases/checkm_db"

python metahict.py reassembly \
  -p "$PWD" --bin out/binning/metahict/metahict_50_10_bins \
  --assembly out/assembly/final_assembly.fasta \
  --hic1 out/preprocessing/final_hic_1.fastq.gz --hic2 out/preprocessing/final_hic_2.fastq.gz \
  --sg1 out/preprocessing/final_sg_1.fastq.gz --sg2 out/preprocessing/final_sg_2.fastq.gz \
  --bam out/alignment/sorted_map.bam --outdir out/reassembly \
  -t 80 -m 128 --checkm2_db "$PWD/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"

python metahict.py scaffolding \
  -p "$PWD" --fasta out/reassembly/reassembled_bins/bin1.fa \
  --hic1 out/preprocessing/final_hic_1.fastq.gz --hic2 out/preprocessing/final_hic_2.fastq.gz \
  --enzyme Sau3AI,MluCI --outdir out/scaffolding -t 80 -m 128 -r 10000 \
  --checkm2_db "$PWD/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"

python metahict.py annotation \
  -p "$PWD" --bin out/reassembly/reassembled_bins --outdir out/annotation -t 80 \
  --gtdbtk_db "$PWD/databases/gtdbtk_db/release220"

python metahict.py mge \
  -p "$PWD" --combined out/reassembly/combined_contigs.fa \
  --contact out/mge_contact/denoised_contact_matrix_normcc.npz \
  --raw-contact out/mge_contact/Raw_contact_matrix.npz \
  --host-taxonomy out/annotation/classify/gtdbtk.bac120.summary.tsv \
  --outdir out/mge -t 80 \
  --genomad_db "$PWD/databases/genomad_db/genomad_db"
```

For exact options, use `python metahict.py <module> --help`.

## Testing

Run the smoke test locally:

```bash
bash nextflow/ci/run_smoke_ci.sh
```

The smoke test runs the DSL2 workflow in Nextflow stub mode and verifies expected output paths with `nextflow/bin/check_expected_outputs.py`. It checks workflow parsing, process connectivity, staging, publication, configuration handling, and expected output structure. 

## Third-party software

METAHICT uses established third-party tools for read processing, assembly, alignment, contact normalization, binning, reassembly support, scaffolding, annotation, and MGE discovery. Sources, versions, licenses, and citations are listed in [THIRD_PARTY.md](THIRD_PARTY.md).

The Python 3-compatible bin3C dependency is maintained as an attributed fork:

```text
https://github.com/1001shiyuan/bin3C-python3/tree/metahict-python3-port
```

## Citation

If you use METAHICT, cite the METAHICT software release and the third-party tools listed in [THIRD_PARTY.md](THIRD_PARTY.md). The METAHICT manuscript describes the benchmark analyses and biological case studies associated with this release.

## License

METAHICT is distributed under the GNU General Public License. See [LICENSE](LICENSE) for details.
