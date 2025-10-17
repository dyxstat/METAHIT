# METAHIT
METAHIT enables comprehensive and flexible genome-resolved microbiome analysis with metagenomic Hi-C.

## Installation
### Dependencies
To install all dependencies required for **METAHIT**, run the setup script located in the `0_installation` folder:

```bash
bash run_setup_in_venv.sh
```

This command will automatically create and activate a minimal Conda environment (`metahit_venv`) and then execute `setup.sh` inside it. During this process, all necessary tools and environments (e.g., BBTools, CheckM2, GTDB-Tk, geNomad, CheckV) will be downloaded, configured, and installed into an `external/` directory within the repository. Once setup completes, you can optionally add `external/bin/` to your system `PATH` for easier access to the installed executables.

### Databases
The folder `0_installation/db_setup` contains four scripts to download and set up databases for **CheckM**, **CheckM2**, **CheckV**, and **GTDB-Tk**. By default, each script downloads the database into a `database/` folder in your current working directory, but you can optionally provide a custom path.

**CheckM database downloading:**  
```bash
bash 0_installation/db_setup/checkm_db.sh [DB_DIR]
```

**CheckM2 database downloading:**  
```bash
bash 0_installation/db_setup/checkm2_db.sh [DB_DIR]
```

**CheckV database downloading:**  
```bash
bash 0_installation/db_setup/checkv_db.sh [DB_DIR]
```

**GTDB-Tk database downloading:**  
```bash
bash 0_installation/db_setup/gtdbtk_db.sh [DB_DIR]
```

## Usage
Once installation and database setup are complete, **METAHIT** can be run by executing each module independently. The framework consists of **10 modules**, each corresponding to a numbered folder in the repository: `1_preprocessing`, `2_assembly`, `3_alignment`, `4_coverage`, `5_contact`, `6_binning`, `7_reassembly`, `8_scaffolding`, `9_annotation`, and `10_MGE`. You can view the overall structure of METAHIT below:

![METAHIT overview](images/Metahit_Overview.png)

### Sample Usage
Each module in **METAHIT** can be executed independently, allowing users to run only specific steps of the workflow. For example, to run the binning module:

```bash
bash 6_binning.sh <FASTA> <BAM> <OUTDIR> <PROJECT_PATH> -t 80
```

**Input arguments:**
- `<FASTA>` — Assembled contigs file  
- `<BAM>` — Mapped Hi-C alignments  
- `<OUTDIR>` — Output directory for binning results  
- `<PROJECT_PATH>` — Path to the METAHIT project directory

Additional optional parameters for fine-tuning MetaCC, bin3C, and ImputeCC are supported. All parameters have internal default values and can be customized as needed; see the source code (`6_binning/6a_binning.py`) for full details.

### Selective Execution
Since the **METAHIT** modules can be executed independently, each step is optional and can be skipped depending on computational resources and analysis needs.

For example, the **reassembly** module is computationally intensive and performs best with sufficient sequencing coverage. In practice, users may choose to reassemble only selected bins—such as those with higher contamination or of particular biological importance—to balance resource use and data quality. When resources are constrained, this step can be skipped, and analyses can proceed using the consolidated bins from the **binning** module, although our benchmarking indicates that reassembly substantially improves contiguity and reduces contamination.

Similarly, the final three modules—**scaffolding**, **annotation**, and **MGE**—are also optional and can be included or omitted depending on the study’s objectives.

## Copyright
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.





























