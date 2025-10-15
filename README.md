# METAHIT
METAHIT enables comprehensive and flexible genome-resolved microbiome analysis with metagenomic Hi-C

## Installation
To install all dependencies required for **METAHIT**, simply run the provided setup script:
```bash
bash run_setup_in_venv.sh
```
This command will automatically create and activate a minimal Conda environment (`metahit_venv`) and then execute `setup.sh` inside it. During this process, all necessary tools and environments (e.g., BBTools, CheckM2, GTDB-Tk, geNomad, CheckV) will be downloaded, configured, and installed into an `external/` directory within the repository. Once setup completes, you can optionally add `external/bin/` to your system `PATH` for easier access to the installed executables.
