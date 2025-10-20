# METAHIT
METAHIT enables comprehensive and flexible genome-resolved microbiome analysis with metagenomic Hi-C.

## Installation
To install **METAHIT**, first clone the repository and navigate to the project directory:

```bash
git clone https://github.com/dyxstat/METAHIT.git
cd METAHIT
```

Once complete, you can proceed to set up dependencies and databases.

### Dependencies
To install all dependencies required for **METAHIT**, run the setup script located in the `installation` folder:

```bash
bash installation/run_setup_in_venv.sh
```

This command will automatically create and activate a minimal Conda environment (`metahit_venv`) and then execute `setup.sh` inside it. During this process, all necessary tools and environments (e.g., BBTools, CheckM, CheckM2, GTDB-Tk, geNomad, CheckV) will be downloaded, configured, and installed into an `external/` directory within the repository. Once setup completes, you can optionally add `external/bin/` to your system `PATH` for easier access to the installed executables.

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
Once installation and database setup are complete, **METAHIT** can be run by executing each module independently. The framework consists of **10 modules**, each corresponding to a numbered folder in the repository: `1_preprocessing`, `2_assembly`, `3_alignment`, `4_coverage`, `5_contact`, `6_binning`, `7_reassembly`, `8_scaffolding`, `9_annotation`, and `10_MGE`. You can view the overall structure of METAHIT below:

![METAHIT overview](images/Metahit_Overview.png)

### Basic Usage
`metahit.py` is the main command-line wrapper that controls all modules in the MetaHit pipeline. Each module can be executed individually using subcommands.

```bash
python metahit.py <module> [options]
```

Example:
```bash
python metahit.py preprocessing -p /path/to/project -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o output_dir -t 80
```

#### 1. Preprocessing Module  
```bash
bash 1_preprocessing.sh -p <PROJECT_PATH> -1 <READ1> -2 <READ2> -o <OUTDIR> [options]
```

**Inputs**  
- <PROJECT_PATH> — Path to the METAHIT project directory  
- <READ1> — Forward shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)  
- <READ2> — Reverse shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)  
- <OUTDIR> — Output directory for preprocessed reads  

**Outputs**  
- <OUTDIR>/final_<prefix>_1.fastq.gz — Final preprocessed forward reads  
- <OUTDIR>/final_<prefix>_2.fastq.gz — Final preprocessed reverse reads  

**Parameters**  
- -t, --threads — Number of CPU threads (default 80)  
- --dedup — Enable duplicate removal for Hi-C reads

#### 2. Assembly Module  
```bash
bash 2_assembly.sh -p <PROJECT_PATH> -1 <READ1> -2 <READ2> -o <OUTDIR> [options]
```

**Inputs**  
- <PROJECT_PATH> — Path to the METAHIT project directory  
- <READ1> — Forward preprocessed shotgun reads (.fastq or .fastq.gz)  
- <READ2> — Reverse preprocessed shotgun reads (.fastq or .fastq.gz)  
- <OUTDIR> — Output directory for assembled contigs  

**Outputs**  
- <OUTDIR>/final_assembly.fasta — Final assembled contigs  

**Parameters**  
- -t, --threads — Number of CPU threads (default 80)  
- -l — Minimum contig length (default 1000 bp)  
- --megahit / --metaspades / --metaflye — Choose assembler (default MEGAHIT)

#### 3. Alignment Module  
```bash
bash 3_alignment.sh -p <PROJECT_PATH> -r <REFERENCE> -1 <READ1> -2 <READ2> -o <OUTDIR> [options]
```

**Inputs**  
- <PROJECT_PATH> — Path to the METAHIT project directory  
- <REFERENCE> — Assembled contigs file (.fasta)  
- <READ1> — Forward preprocessed Hi-C reads (.fastq or .fastq.gz)  
- <READ2> — Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)  
- <OUTDIR> — Output directory for alignment results  

**Outputs**  
- <OUTDIR>/sorted_map.bam — Sorted BAM file of aligned Hi-C reads  
- <OUTDIR>/3d_ratio.txt — 3D ratio 

**Parameters**  
- -t, --threads — Number of CPU threads (default 80)  
- --samtools-filter — Filtering flag for `samtools view` (default `-F 0x900`)

#### 4. Coverage Module  
```bash
bash 4_coverage.sh -p <PROJECT_PATH> -1 <READ1> -2 <READ2> -r <REFERENCE> -o <OUTDIR> [options]
```

**Inputs**  
- <PROJECT_PATH> — Path to the METAHIT project directory  
- <READ1> — Forward shotgun reads (.fastq or .fastq.gz)  
- <READ2> — Reverse shotgun reads (.fastq or .fastq.gz)  
- <REFERENCE> — Assembled contigs file (.fasta)  
- <OUTDIR> — Output directory for coverage results  

**Outputs**  
- <OUTDIR>/SG_map_sorted.bam — Sorted BAM file of mapped reads  
- <OUTDIR>/coverage.txt — Contig-level coverage summary  
- <OUTDIR>/pair.txt — Paired-contig information  

**Parameters**  
- -t, --threads — Number of CPU threads (default 80)

#### 5. Contact Module  
```bash
bash 5_contact.sh <METHOD> -p <PROJECT_PATH> --bam <BAM> --fasta <FASTA> --out <OUTDIR> --enzyme <ENZYME>
```

**Inputs**  
- <METHOD> — Normalization method (`metator`, `hiczin`, `normcc`, etc.)  
- <PROJECT_PATH> — Path to the METAHIT project directory  
- <BAM> — Hi-C read alignment file (.bam)  
- <FASTA> — Assembled contigs file (.fasta)  
- <OUTDIR> — Output directory for contact maps  
- <ENZYME> — Restriction enzymes used in the Hi-C library (e.g., Sau3AI,MluCI)  

**Outputs**  
- <OUTDIR>/Raw_contact_matrix.npz — Raw contact matrix   
- <OUTDIR>/normalized_contact_matrix.npz — Normalized contact matrix
- <OUTDIR>/contig_info.csv — Contig metadata 

#### 6. Binning Module  
```bash
bash 6_binning.sh <FASTA> <BAM> <OUTDIR> <PROJECT_PATH> [options]  
```

**Inputs**  
- `<FASTA>` — Assembled contigs file (.fa or .fasta)  
- `<BAM>` — Hi-C reads aligned to the contigs (.bam) 
- `<OUTDIR>` — Output directory for binning results 
- `<PROJECT_PATH>` — Path to the METAHIT project directory 

**Outputs**  
- `<OUTDIR>/bin3c/fasta/*.fna` — bin3C bins 
- `<OUTDIR>/metacc/BIN/*.fa` — MetaCC bins 
- `<OUTDIR>/imputecc/FINAL_BIN/*.fa` — ImputeCC bins 
- `<OUTDIR>/metahit/metahit_50_10_bins/*.fa` — Integrated final bins produced by MetaHIT bin refinement  
- `<OUTDIR>/metahit/figures/heatmap.png` — Heatmap of final bins  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)
- `--checkm_db` - Custom path for the CheckM database

#### 7. Reassembly Module
```bash
bash 7_reassembly.sh -p <PROJECT_PATH> --bin <BIN_DIR> --assembly <ASSEMBLY> --hic1 <HIC_READ1> --hic2 <HIC_READ2> --sg1 <SHOTGUN_READ1> --sg2 <SHOTGUN_READ2> --bam <BAM> --outdir <OUTDIR> -t <THREADS> -m <MEMORY>
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHIT project directory  
- `<BIN_DIR>` — Directory containing input bins  
- `<ASSEMBLY>` — Original assembly FASTA file (.fa or .fasta)  
- <HIC_READ1> — Forward preprocessed Hi-C reads (.fastq or .fastq.gz)  
- <HIC_READ2> — Reverse preprocessed Hi-C reads (.fastq or .fastq.gz) 
- <SHOTGUN_READ1> — Forward preprocessed shotgun reads (.fastq or .fastq.gz)  
- <SHOTGUN_READ2> — Reverse preprocessed shotgun reads (.fastq or .fastq.gz)  
- `<BAM>` — Hi-C read alignments to the assembly (.bam)  
- `<OUTDIR>` — Output directory for reassembly results  

**Outputs**  
- `<OUTDIR>/reassembled_bins/` — Final reassembled bins (.fa)  
- `<OUTDIR>/unmapped_assembly/final.contigs.fa` — Assembly of unmapped reads  
- `<OUTDIR>/combined/combined_contigs.fa` — Combined contigs (bins and unmapped)  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `-m, --memory` — Memory in GB (default 24)  
- `--parallel` — Enable per-bin parallel reassembly (1 thread per bin)
- `--checkm2_db` - Custom path for the CheckM2 database

#### 8. Scaffolding Module  
```bash
bash 8_scaffolding.sh -p <PROJECT_PATH> --fasta <BIN_FASTA> --bam <BAM> --enzyme <ENZYME> --outdir <OUTDIR> --hic1 <HIC1> --hic2 <HIC2> -t <THREADS> -m <MEMORY> -r <RESOLUTION>
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHIT project directory  
- `<BIN_FASTA>` — Input bin file for scaffolding (.fa or .fasta)  
- `<BAM>` — Optional Hi-C read alignments to the assembly (.bam) 
- `<ENZYME>` — Restriction enzymes used in Hi-C library (e.g., Sau3AI,MluCI)  
- `<OUTDIR>` — Output directory for scaffolding results  
- <HIC1> — Forward preprocessed Hi-C reads (.fastq or .fastq.gz)  
- <HIC2> — Reverse preprocessed Hi-C reads (.fastq or .fastq.gz) 

**Outputs**  
- `<OUTDIR>/yahs/scaffold_scaffolds_final.fa` — Scaffolded genome from YaHS  
- `<OUTDIR>/figures/heatmap.png` — Contact heatmap of scaffolded genome  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)  
- `-m, --memory` — Memory limit for YaHS and SPAdes (default: 80% of available RAM)  
- `-r, --resolution` — Segment length for visualization (default 1000 bp)  
- `--bam` — Skip new Hi-C alignment by providing an existing BAM
- `--checkm2_db` - Custom path for the CheckM2 database

#### 9. Annotation Module
```bash
bash 9_annotation.sh -p <PROJECT_PATH> --bin <BIN_DIR> --outdir <OUTDIR> -t <THREADS>  
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHIT project directory  
- `<BIN_DIR>` — Directory containing input bins 
- `<OUTDIR>` — Output directory for annotation results  

**Outputs**  
- `<OUTDIR>/gtdbtk.bac120.summary.tsv` — Summary of GTDB-Tk bacterial classifications  
- `<OUTDIR>/gtdbtk.ar122.summary.tsv` — Summary of GTDB-Tk archaeal classifications  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)
- `--gtdbtk_db` - Custom path for the GTDB-Tk database

#### 10. MGE Module  
```bash
bash 10_MGE.sh -p <PROJECT_PATH> --combined <COMBINED_FASTA> --contact <CONTACT_MATRIX> --outdir <OUTDIR> -t <THREADS>  
```

**Inputs**  
- `<PROJECT_PATH>` — Path to the METAHIT project directory  
- `<COMBINED_FASTA>` — Combined contigs FASTA file (include both binned and unmapped contigs, .fa)  
- `<CONTACT_MATRIX>` — Normalized contact matrix (.npz) 
- `<OUTDIR>` — Output directory for MGE analysis results  

**Outputs**  
- `<OUTDIR>/genomad_output/` — geNomad predictions of viral and plasmid contigs  
- `<OUTDIR>/checkv_output/virus/quality_summary.tsv` — CheckV QC summary of viral contigs  

**Parameters**  
- `-t, --threads` — Number of CPU threads (default 80)
- `--genomad_db` - Custom path for the geNomad database
- `--checkv_db` - Custom path for the CheckV database

#### Selective Execution
Since the **METAHIT** modules can be executed independently, each step is optional and can be skipped depending on computational resources and analysis needs.

For example, the **reassembly** module is computationally intensive and performs best with sufficient sequencing coverage. In practice, users may choose to reassemble only selected bins—such as those with higher contamination or of particular biological importance—to balance resource use and data quality. When resources are constrained, this step can be skipped, and analyses can proceed using the consolidated bins from the **binning** module, although our benchmarking indicates that reassembly substantially improves contiguity and reduces contamination.

Similarly, the final three modules—**scaffolding**, **annotation**, and **MGE**—are also optional and can be included or omitted depending on the study’s objectives.

## Copyright
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.





















