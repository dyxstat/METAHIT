#!/usr/bin/env python
import argparse
import subprocess
import os
import sys

script_dir = sys.path[0]

def run_command(command):
    print(f"[INFO] Executing command:\n{command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {e}")
        sys.exit(1)

def ensure_dir_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def absolute_path(path):
    return os.path.abspath(path)

def preprocessing(args):
    print("[INFO] Running Preprocessing Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/1_preprocessing/1_preprocessing.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )
    
    if args.prefix:
        command += f' --prefix "{args.prefix}"'
    if args.dedup:
        command += " --dedup"

    run_command(command)
    
def assembly(args):
    print("[INFO] Running Assembly Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/2_assembly/2_assembly.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.min_len:
        command += f' -l {args.min_len}'
    if args.megahit:
        command += " --megahit"
    if args.metaspades:
        command += " --metaspades"
    if args.metaflye:
        command += " --metaflye"

    run_command(command)
    
def alignment(args):
    print("[INFO] Running Alignment Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/3_alignment/3_alignment.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-r "{absolute_path(args.reference)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.samtools_filter:
        command += f' --samtools-filter "{args.samtools_filter}"'

    run_command(command)
    
def coverage(args):
    print("[INFO] Running Coverage Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/4_coverage/4_coverage.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-r "{absolute_path(args.reference)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )

    run_command(command)
    
def contact(args):
    print("[INFO] Running Contact Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/5_contact/5_contact.sh" '
        f'{args.method} '
        f'-p "{absolute_path(args.project_path)}" '
        f'--bam "{absolute_path(args.bam)}" '
        f'--fasta "{absolute_path(args.fasta)}" '
        f'--out "{output_dir}" '
        f'--enzyme "{args.enzyme}"'
    )

    run_command(command)
    
def binning(args):
    print("[INFO] Running Binning Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/6_binning/6_binning.sh" '
        f'"{absolute_path(args.fasta)}" '
        f'"{absolute_path(args.bam)}" '
        f'"{output_dir}" '
        f'"{absolute_path(args.project_path)}" '
        f'--threads {args.threads}'
    )
    
    if args.checkm_db:
        command += f' --checkm_db "{absolute_path(args.checkm_db)}"'

    run_command(command)
    
def reassembly(args):
    print("[INFO] Running Reassembly Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/7_reassembly/7_reassembly.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--bin "{absolute_path(args.bin)}" '
        f'--assembly "{absolute_path(args.assembly)}" '
        f'--hic1 "{absolute_path(args.hic1)}" '
        f'--hic2 "{absolute_path(args.hic2)}" '
        f'--sg1 "{absolute_path(args.sg1)}" '
        f'--sg2 "{absolute_path(args.sg2)}" '
        f'--bam "{absolute_path(args.bam)}" '
        f'--outdir "{output_dir}" '
        f'-t {args.threads} '
        f'-m {args.memory}'
    )

    if args.parallel:
        command += " --parallel"
    if args.checkm2_db:
        command += f' --checkm2_db "{absolute_path(args.checkm2_db)}"'

    run_command(command)
    
def scaffolding(args):
    print("[INFO] Running Scaffolding Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/8_scaffolding/8_scaffolding.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--fasta "{absolute_path(args.fasta)}" '
        f'--enzyme "{args.enzyme}" '
        f'--outdir "{output_dir}" '
        f'--hic1 "{absolute_path(args.hic1)}" '
        f'--hic2 "{absolute_path(args.hic2)}" '
        f'-t {args.threads} '
        f'-r {args.resolution}'
    )

    if args.bam:
        command += f' --bam "{absolute_path(args.bam)}"'
    if args.memory:
        command += f' -m {args.memory}'    
    if args.checkm2_db:
        command += f' --checkm2_db "{absolute_path(args.checkm2_db)}"'
        
    run_command(command)
    
def annotation(args):
    print("[INFO] Running Annotation Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/9_annotation/9_annotation.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--bin "{absolute_path(args.bin)}" '
        f'--outdir "{output_dir}" '
        f'-t {args.threads}'
    )
    
    if args.gtdbtk_db:
        command += f' --gtdbtk_db "{absolute_path(args.gtdbtk_db)}"'

    run_command(command)

def mge(args):
    print("[INFO] Running MGE Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/10_MGE/10_MGE.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--combined "{absolute_path(args.combined)}" '
        f'--contact "{absolute_path(args.contact)}" '
        f'--outdir "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.genomad_db:
        command += f' --genomad_db "{absolute_path(args.genomad_db)}"'
    if args.checkv_db:
        command += f' --checkv_db "{absolute_path(args.checkv_db)}"'

    run_command(command)

def main():
    parser = argparse.ArgumentParser(description="MetaHit Pipeline Wrapper")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Preprocessing subcommand 
    pre_parser = subparsers.add_parser("preprocessing", help="Run preprocessing of shotgun or Hi-C reads")
    pre_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    pre_parser.add_argument("-1", dest="r1", required=True, help="Forward shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)")
    pre_parser.add_argument("-2", dest="r2", required=True, help="Reverse shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)")
    pre_parser.add_argument("-o", "--output", required=True, help="Output directory for preprocessed reads")
    pre_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    pre_parser.add_argument("--dedup", action="store_true", help="Enable duplicate removal for Hi-C reads")
    pre_parser.add_argument("--prefix", help="Custom prefix for output files (default=derived from input filename)")
    pre_parser.set_defaults(func=preprocessing)

    # Assembly subcommand
    asm_parser = subparsers.add_parser("assembly", help="Run assembly module to generate contigs")
    asm_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    asm_parser.add_argument("-1", dest="r1", required=True, help="Forward preprocessed reads (.fastq or .fastq.gz)")
    asm_parser.add_argument("-2", dest="r2", required=True, help="Reverse preprocessed reads (.fastq or .fastq.gz)")
    asm_parser.add_argument("-o", "--output", required=True, help="Output directory for assembled contigs")
    asm_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    asm_parser.add_argument("-l", "--min-len", type=int, default=1000, help="Minimum contig length (default=1000 bp)")
    asm_parser.add_argument("--megahit", action="store_true", help="Assemble with MEGAHIT (default)")
    asm_parser.add_argument("--metaspades", action="store_true", help="Assemble with metaSPAdes")
    asm_parser.add_argument("--metaflye", action="store_true", help="Assemble with metaFlye")
    asm_parser.set_defaults(func=assembly)

    # Alignment subcommand
    aln_parser = subparsers.add_parser("alignment", help="Run alignment module for Hi-C read mapping")
    aln_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    aln_parser.add_argument("-r", "--reference", required=True, help="Assembled contigs file (.fasta)")
    aln_parser.add_argument("-1", dest="r1", required=True, help="Forward preprocessed Hi-C reads (.fastq or .fastq.gz)")
    aln_parser.add_argument("-2", dest="r2", required=True, help="Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)")
    aln_parser.add_argument("-o", "--output", required=True, help="Output directory for alignment results")
    aln_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    aln_parser.add_argument("--samtools-filter", default="-F 0x900", help="Filtering flag for samtools view (default='-F 0x900')")
    aln_parser.set_defaults(func=alignment)
    
    # Coverage subcommand
    cov_parser = subparsers.add_parser("coverage", help="Run coverage module for read mapping and coverage estimation")
    cov_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    cov_parser.add_argument("-1", dest="r1", required=True, help="Forward shotgun reads (.fastq or .fastq.gz)")
    cov_parser.add_argument("-2", dest="r2", required=True, help="Reverse shotgun reads (.fastq or .fastq.gz)")
    cov_parser.add_argument("-r", "--reference", required=True, help="Assembled contigs file (.fasta)")
    cov_parser.add_argument("-o", "--output", required=True, help="Output directory for coverage results")
    cov_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    cov_parser.set_defaults(func=coverage)
    
    # Contact subcommand
    con_parser = subparsers.add_parser("contact", help="Run contact module to generate contact matrices")
    con_parser.add_argument("method", help="Normalization method (e.g., metator, hiczin, normcc)")
    con_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    con_parser.add_argument("--bam", required=True, help="Hi-C read alignment file (.bam)")
    con_parser.add_argument("--fasta", required=True, help="Assembled contigs file (.fasta)")
    con_parser.add_argument("--out", required=True, help="Output directory for contact maps")
    con_parser.add_argument("--enzyme", required=True, help="Restriction enzymes used in Hi-C library (e.g., Sau3AI,MluCI)")
    con_parser.set_defaults(func=contact)
    
    # Binning subcommand
    bin_parser = subparsers.add_parser("binning", help="Run binning module to generate MAGs and refined bins")
    bin_parser.add_argument("--fasta", required=True, help="Assembled contigs file (.fa or .fasta)")
    bin_parser.add_argument("--bam", required=True, help="Hi-C reads aligned to the contigs (.bam)")
    bin_parser.add_argument("--output", required=True, help="Output directory for binning results")
    bin_parser.add_argument("--project_path", required=True, help="Path to the MetaHit project directory")
    bin_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    bin_parser.add_argument("--checkm_db", help="Path to CheckM database (if not using default environment variable)")
    bin_parser.set_defaults(func=binning)
    
    # Reassembly subcommand
    reas_parser = subparsers.add_parser("reassembly", help="Run reassembly module to reconstruct bins and unmapped reads")
    reas_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    reas_parser.add_argument("--bin", required=True, help="Directory containing input bins")
    reas_parser.add_argument("--assembly", required=True, help="Original assembly FASTA file (.fa or .fasta)")
    reas_parser.add_argument("--hic1", required=True, help="Forward preprocessed Hi-C reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--hic2", required=True, help="Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--sg1", required=True, help="Forward preprocessed shotgun reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--sg2", required=True, help="Reverse preprocessed shotgun reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--bam", required=True, help="Hi-C read alignments to the assembly (.bam)")
    reas_parser.add_argument("--outdir", required=True, help="Output directory for reassembly results")
    reas_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default 80)")
    reas_parser.add_argument("-m", "--memory", type=int, default=24, help="Memory in GB (default 24)")
    reas_parser.add_argument("--parallel", action="store_true", help="Enable per-bin parallel reassembly (1 thread per bin)")
    reas_parser.add_argument("--checkm2_db", help="Path to CheckM2 database (if not using default environment variable)")
    reas_parser.set_defaults(func=reassembly)

    # Scaffolding subcommand
    scaf_parser = subparsers.add_parser("scaffolding", help="Run scaffolding module to generate scaffolded genomes and heatmaps")
    scaf_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    scaf_parser.add_argument("--fasta", required=True, help="Input bin FASTA file for scaffolding (.fa or .fasta)")
    scaf_parser.add_argument("--bam", help="Optional Hi-C read alignments to the assembly (.bam)")
    scaf_parser.add_argument("--enzyme", required=True, help="Restriction enzymes used in the Hi-C library (e.g., Sau3AI,MluCI)")
    scaf_parser.add_argument("--outdir", required=True, help="Output directory for scaffolding results")
    scaf_parser.add_argument("--hic1", required=True, help="Forward preprocessed Hi-C reads (.fastq or .fastq.gz)")
    scaf_parser.add_argument("--hic2", required=True, help="Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)")
    scaf_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    scaf_parser.add_argument("-m", "--memory", type=str, default=None, help="Memory limit (default=80% of available RAM)")
    scaf_parser.add_argument("-r", "--resolution", type=int, default=10000, help="Segment length for visualization (default=10000 bp)")
    scaf_parser.add_argument("--checkm2_db", help="Path to CheckM2 database (if not using default environment variable)")
    scaf_parser.set_defaults(func=scaffolding)

    # Annotation subcommand
    anno_parser = subparsers.add_parser("annotation", help="Run annotation module using GTDB-Tk for taxonomic classification")
    anno_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    anno_parser.add_argument("--bin", required=True, help="Directory containing input bins")
    anno_parser.add_argument("--outdir", required=True, help="Output directory for annotation results")
    anno_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    anno_parser.add_argument("--gtdbtk_db", help="Path to GTDB-Tk database (if not using default environment variable)")
    anno_parser.set_defaults(func=annotation)

    # MGE subcommand
    mge_parser = subparsers.add_parser("mge", help="Run MGE module for viral/plasmid detection and host–MGE interaction analysis")
    mge_parser.add_argument("-p", "--project_path", required=True, help="Path to the MetaHit project directory")
    mge_parser.add_argument("--combined", required=True, help="Combined contigs FASTA file (includes binned and unmapped contigs)")
    mge_parser.add_argument("--contact", required=True, help="Normalized contact matrix (.npz)")
    mge_parser.add_argument("--outdir", required=True, help="Output directory for MGE analysis results")
    mge_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    mge_parser.add_argument("--genomad_db", help="Path to geNomad database (if not using default environment variable)")
    mge_parser.add_argument("--checkv_db", help="Path to CheckV database (if not using default environment variable)")
    mge_parser.set_defaults(func=mge)

    args = parser.parse_args()
    if args.command == "preprocessing":
        preprocessing(args)
        
    elif args.command == "assembly":
        assembly(args)
        
    elif args.command == "alignment":
        alignment(args)

    elif args.command == "coverage":
        coverage(args)
        
    elif args.command == "contact":
        contact(args)
        
    elif args.command == "binning":
        binning(args)
        
    elif args.command == "reassembly":
        reassembly(args)

    elif args.command == "scaffolding":
        scaffolding(args)

    elif args.command == "annotation":
        annotation(args)

    elif args.command == "mge":
        mge(args)

if __name__ == "__main__":
    main()
