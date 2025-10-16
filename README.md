# METAHIT
METAHIT enables comprehensive and flexible genome-resolved microbiome analysis with metagenomic Hi-C

## Installation
### Dependencies
To install all dependencies required for **METAHIT**, run the setup script located in the `0_installation` folder:
```bash
bash run_setup_in_venv.sh
```
This command will automatically create and activate a minimal Conda environment (`metahit_venv`) and then execute `setup.sh` inside it. During this process, all necessary tools and environments (e.g., BBTools, CheckM2, GTDB-Tk, geNomad, CheckV) will be downloaded, configured, and installed into an `external/` directory within the repository. Once setup completes, you can optionally add `external/bin/` to your system `PATH` for easier access to the installed executables.

### Databases
The folder `0_installation/db_setup` contains four scripts to download and set up databases for **CheckM**, **CheckM2**, **CheckV**, and **GTDB-Tk**.  
By default, each script downloads the database into a `database/` folder in your current working directory, but you can optionally provide a custom path.
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
