# METAHICT
METAHICT enables comprehensive and flexible genome-resolved microbiome analysis with metagenomic Hi-C.

## Installation
To install **METAHICT**, first clone the repository and navigate to the project directory:

```bash
git clone https://github.com/dyxstat/METAHICT.git
cd METAHICT
```

Once complete, you can proceed to set up dependencies and databases.

### Dependencies
To install all dependencies required for **METAHICT**, run the setup script located in the `installation` folder:

```bash
bash installation/run_setup_in_venv.sh
```

This command creates a minimal bootstrap Conda environment (`metahict_venv`) and
then executes `setup.sh`.  The scientific software is installed into the
project-local `conda_envs/` environments from the exact Linux lock files.
BBTools is the exception: it is downloaded to `external/bbmap/`, because the
preprocessing and coverage scripts call its shell entry points directly.  Once
setup completes, you can optionally add `external/bin/` to your system `PATH`.

The Linux release installs each Conda environment from an exact package-artifact
lock under `installation/locks/linux-64/`; it does not re-solve floating
dependencies. `installation/env.yaml` and `installation/checkm2.yaml` are the
human-readable specifications, while the lock files are the authoritative
reproducible installation inputs.

The release configuration installs the Python 3-compatible bin3C dependency
from the separately maintained, AGPL-3.0-licensed
`1001shiyuan/bin3C-python3` fork at the immutable commit pinned in
`installation/env.yaml`. Its upstream provenance, licence, and citation are
recorded in [THIRD_PARTY.md](THIRD_PARTY.md).

### Nextflow workflow

The native DSL2 Nextflow workflow is in [`nextflow/`](nextflow/).  It calls the
same numbered METAHICT modules described below, while Nextflow controls their
order, inputs, outputs, resources, provenance reports, and restart behavior.
Users choose one dependency-execution profile: project-local Conda, Docker, or
Apptainer.  Docker and Apptainer use a separately published container image
once the tested release image is available; the image is not stored in this
Git repository.  See [nextflow/README.md](nextflow/README.md) and
[nextflow/REPRODUCIBILITY.md](nextflow/REPRODUCIBILITY.md) for exact workflow
behavior, defaults, environments, and database requirements.

### Databases
The folder `installation/db` contains five scripts to download and set up databases for **CheckM**, **CheckM2**, **GTDB-Tk**, **geNomad** and **CheckV**. By default, each script downloads the database into a `databases/` folder in your current working directory. You can optionally specify a custom path during installation. If you do, please make sure to provide the same custom path when running the corresponding modules — CheckM for module 6, CheckM2 for modules 7 and 8, GTDB-Tk for module 9, and geNomad and CheckV for module 10.

**CheckM database downloading:**  
```bash
bash installation/db/checkm_db.sh [DB_DIR]
```

**CheckM2 database downloading:**  
```bash
bash installation/db/checkm2_db.sh [DB_DIR]
```

**GTDB-Tk database downloading:**  
```bash
bash installation/db/gtdbtk_db.sh [DB_DIR]
```

**geNomad database downloading:**  
```bash
bash installation/db/genomad_db.sh [DB_DIR]
```

**CheckV database downloading:**  
```bash
bash installation/db/checkv_db.sh [DB_DIR]
```

## Usage
Once installation and database setup are complete, **METAHICT** can be run by executing each module independently. The framework consists of **10 modules**, each corresponding to a numbered folder in the repository: `1_preprocessing`, `2_assembly`, `3_alignment`, `4_coverage`, `5_contact`, `6_binning`, `7_reassembly`, `8_scaffolding`, `9_annotation`, and `10_MGE`. You can view the overall structure of METAHICT below:

![METAHICT overview](images/METAHICT_Overview.png)

### Basic Usage
`metahict.py` is the main command-line wrapper that controls all modules in the METAHICT pipeline. Each module can be executed individually using subcommands.

```bash
conda activate metahict_env
python metahict.py <module> [options]
```

### 1. Preprocessing Module  
```bash
python metahict.py preprocessing -p <PROJECT_PATH> -1 <READ1> -2 <READ2> -o <OUTDIR> [options]
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<READ1>` — Forward shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)  
- `<READ2>` — Reverse shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)  
- `<OUTDIR>` — Output directory for preprocessed reads  

**Outputs**  
- `<OUTDIR>/final_<prefix>_1.fastq.gz` — Final preprocessed forward reads  
- `<OUTDIR>/final_<prefix>_2.fastq.gz` — Final preprocessed reverse reads  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `--dedup` — Enable duplicate removal for Hi-C reads
- `--prefix` — Custom prefix for output files (default: base name derived from input reads)

### 2. Assembly Module  
```bash
python metahict.py assembly -p <PROJECT_PATH> -1 <READ1> -2 <READ2> -o <OUTDIR> [options]
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<READ1>` — Forward preprocessed shotgun reads (.fastq or .fastq.gz)  
- `<READ2>` — Reverse preprocessed shotgun reads (.fastq or .fastq.gz)  
- `<OUTDIR>` — Output directory for assembled contigs  

**Outputs**  
- `<OUTDIR>/final_assembly.fasta` — Final assembled contigs  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `-l` — Minimum contig length (default 1000 bp)  
- `--megahit / --metaspades / --metaflye` — Choose assembler (default MEGAHIT)

### 3. Alignment Module  
```bash
python metahict.py alignment -p <PROJECT_PATH> -r <REFERENCE> -1 <READ1> -2 <READ2> -o <OUTDIR> [options]
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<REFERENCE>` — Assembled contigs file (.fasta)  
- `<READ1>` — Forward preprocessed Hi-C reads (.fastq or .fastq.gz)  
- `<READ2>` — Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)  
- `<OUTDIR>` — Output directory for alignment results  

**Outputs**  
- `<OUTDIR>/sorted_map.bam` — Sorted BAM file of aligned Hi-C reads  
- `<OUTDIR>/3d_ratio.txt` — 3D ratio 

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `--samtools-filter` — Filtering flag for `samtools view` (default `-F 0x900`)

### 4. Coverage Module  
```bash
python metahict.py coverage -p <PROJECT_PATH> -1 <READ1> -2 <READ2> -r <REFERENCE> -o <OUTDIR> [options]
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<READ1>` — Forward shotgun reads (.fastq or .fastq.gz)  
- `<READ2>` — Reverse shotgun reads (.fastq or .fastq.gz)  
- `<REFERENCE>` — Assembled contigs file (.fasta)  
- `<OUTDIR>` — Output directory for coverage results  

**Outputs**  
- `<OUTDIR>/SG_map_sorted.bam` — Sorted BAM file of mapped reads  
- `<OUTDIR>/coverage.txt` — Contig-level coverage summary  
- `<OUTDIR>/pair.txt` — Paired-contig information  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)

### 5. Contact Module  
```bash
python metahict.py contact <METHOD> -p <PROJECT_PATH> --bam <BAM> --fasta <FASTA> --out <OUTDIR> --enzyme <ENZYME>
```

**Inputs**  
- `<METHOD>` — Normalization method (`metator`, `hiczin`, `normcc`, etc.)  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<BAM>` — Hi-C read alignment file (.bam)  
- `<FASTA>` — Assembled contigs file (.fasta)  
- `<OUTDIR>` — Output directory for contact maps  
- `<ENZYME>` — Restriction enzymes used in the Hi-C library (e.g., Sau3AI,MluCI)  

**Outputs**  
- `<OUTDIR>/Raw_contact_matrix.npz` — Raw contact matrix   
- `<OUTDIR>/normalized_contact_matrix.npz` — Normalized contact matrix
- `<OUTDIR>/contig_info.csv` — Contig metadata 

### 6. Binning Module  
```bash
python metahict.py binning <FASTA> <BAM> <OUTDIR> <PROJECT_PATH> [options]  
```

**Inputs**  
- `<FASTA>` — Assembled contigs file (.fa or .fasta)  
- `<BAM>` — Hi-C reads aligned to the contigs (.bam) 
- `<OUTDIR>` — Output directory for binning results 
- `<PROJECT_PATH>` — Path to the METAHICT project directory 

**Outputs**  
- `<OUTDIR>/bin3c/fasta/*.fna` — bin3C bins 
- `<OUTDIR>/metacc/BIN/*.fa` — MetaCC bins 
- `<OUTDIR>/imputecc/FINAL_BIN/*.fa` — ImputeCC bins 
- `<OUTDIR>/metahict/metahict_50_10_bins/*.fa` — Integrated final bins produced by MetaHIT bin refinement  
- `<OUTDIR>/metahict/figures/heatmap.png` — Heatmap of final bins  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)
- `--checkm_db` - Custom path for the CheckM database

### 7. Reassembly Module
```bash
python metahict.py reassembly -p <PROJECT_PATH> --bin <BIN_DIR> --assembly <ASSEMBLY> --hic1 <HIC_READ1> --hic2 <HIC_READ2> --sg1 <SHOTGUN_READ1> --sg2 <SHOTGUN_READ2> --bam <BAM> --outdir <OUTDIR> -t <THREADS> -m <MEMORY>
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<BIN_DIR>` — Directory containing input bins  
- `<ASSEMBLY>` — Original assembly FASTA file (.fa or .fasta)  
- `<HIC_READ1>` — Forward preprocessed Hi-C reads (.fastq or .fastq.gz)  
- `<HIC_READ2>` — Reverse preprocessed Hi-C reads (.fastq or .fastq.gz) 
- `<SHOTGUN_READ1>` — Forward preprocessed shotgun reads (.fastq or .fastq.gz)  
- `<SHOTGUN_READ2>` — Reverse preprocessed shotgun reads (.fastq or .fastq.gz)  
- `<BAM>` — Hi-C read alignments to the assembly (.bam)  
- `<OUTDIR>` — Output directory for reassembly results  

**Outputs**  
- `<OUTDIR>/reassembled_bins/` — Final selected bins as `.fa` files named with clean bin IDs, for example `bin113.fa`  
- `<OUTDIR>/reassembled_bin_name_map.tsv` — Mapping from each clean final bin ID to the selected `.orig`, `.strict`, or `.permissive` candidate  
- `<OUTDIR>/unmapped_assembly/final.contigs.fa` — Assembly of unmapped reads  
- `<OUTDIR>/combined/combined_contigs.fa` — Combined contigs (bins and unmapped)  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `-m, --memory` — Memory in GB (default 24)  
- `--parallel` — Enable per-bin parallel reassembly (1 thread per bin)
- `--checkm2_db` - Custom path for the CheckM2 database

### 8. Scaffolding Module  
```bash
python metahict.py scaffolding -p <PROJECT_PATH> --fasta <BIN_FASTA> --bam <BAM> --enzyme <ENZYME> --outdir <OUTDIR> --hic1 <HIC1> --hic2 <HIC2> -t <THREADS> -m <MEMORY> -r <RESOLUTION>
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<BIN_FASTA>` — Input bin file for scaffolding (.fa or .fasta)  
- `<BAM>` — Optional Hi-C read alignments to the assembly (.bam) 
- `<ENZYME>` — Restriction enzymes used in Hi-C library (e.g., Sau3AI,MluCI)  
- `<OUTDIR>` — Output directory for scaffolding results  
- `<HIC1>` — Forward preprocessed Hi-C reads (.fastq or .fastq.gz)  
- `<HIC2>` — Reverse preprocessed Hi-C reads (.fastq or .fastq.gz) 

**Outputs**  
- `<OUTDIR>/yahs/scaffold_scaffolds_final.fa` — Scaffolded genome from YaHS  
- `<OUTDIR>/figures/heatmap.png` — Contact heatmap of scaffolded genome  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `-m, --memory` — Memory limit for YaHS and SPAdes (default: 80% of available RAM)  
- `-r, --resolution` — Segment length for visualization (default 1000 bp)  
- `--bam` — Skip new Hi-C alignment by providing an existing BAM
- `--checkm2_db` - Custom path for the CheckM2 database

### 9. Annotation Module
```bash
python metahict.py annotation -p <PROJECT_PATH> --bin <BIN_DIR> --outdir <OUTDIR> -t <THREADS>  
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<BIN_DIR>` — Directory containing input bins 
- `<OUTDIR>` — Output directory for annotation results  

**Outputs**  
- `<OUTDIR>/gtdbtk.bac120.summary.tsv` — Summary of GTDB-Tk bacterial classifications  
- `<OUTDIR>/gtdbtk.ar122.summary.tsv` — Summary of GTDB-Tk archaeal classifications  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)
- `--gtdbtk_db` - Custom path for the GTDB-Tk database

### 10. MGE Module  
```bash
python metahict.py mge -p <PROJECT_PATH> --combined <COMBINED_FASTA> --contact <CONTACT_MATRIX> --raw-contact <RAW_CONTACT_MATRIX> --outdir <OUTDIR> -t <THREADS>
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHICT project directory  
- `<COMBINED_FASTA>` — Combined contigs FASTA file (include both binned and unmapped contigs, .fa)  
- `<CONTACT_MATRIX>` — Normalized contact matrix (.npz)
- `<RAW_CONTACT_MATRIX>` — Raw Hi-C contact matrix (.npz), used to count raw read-pair support
- `<OUTDIR>` — Output directory for MGE analysis results  

**Outputs**  
- `<OUTDIR>/candidate_mge_mag_associations_zscore_filtered.tsv` — Filtered candidate viral MGE-MAG association table
- `<OUTDIR>/genomad_output/combined_contigs_summary/combined_contigs_virus_summary.tsv` — geNomad summary of viral contigs
- `<OUTDIR>/genomad_output/combined_contigs_summary/combined_contigs_plasmid_summary.tsv` — geNomad summary of plasmid contigs
- `<OUTDIR>/checkv_output/virus/quality_summary.tsv` — CheckV QC summary of viral contigs  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)
- `--genomad_db` - Custom path for the geNomad database
- `--checkv_db` - Custom path for the CheckV database

### Selective Execution
Since the **METAHICT** modules can be executed independently, each step is optional and can be skipped depending on computational resources and analysis needs.

For example, the **reassembly** module is computationally intensive and performs best with sufficient sequencing coverage. In practice, users may choose to reassemble only selected bins—such as those with higher contamination or of particular biological importance—to balance resource use and data quality. When resources are constrained, this step can be skipped, and analyses can proceed using the consolidated bins from the **binning** module, although our benchmarking indicates that reassembly substantially improves contiguity and reduces contamination.

Similarly, the final three modules—**scaffolding**, **annotation**, and **MGE**—are also optional and can be included or omitted depending on the study’s objectives.

## Copyright
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.















