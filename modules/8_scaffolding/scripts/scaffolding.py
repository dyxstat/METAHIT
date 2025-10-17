#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scaffolding Pipeline for Metagenomic Genome Assembly

This script implements a comprehensive scaffolding workflow for metagenomic bin refinement:
- Filters bin FASTA to remove short contigs
- Runs YaHS for genomic scaffolding (supports both input BAM or fresh alignment)
- Cuts scaffolds into fixed-length segments for visualization
- Generates and normalizes contact matrices
- Performs Hi-C read realignment to scaffolded genome
- Creates scaffold mapping for visualization
- Generates contact heatmap
- Calculates scaffolding quality metrics (N50/L50, CheckM2, Hi-C enrichment)

Pipeline Steps:
1. Filter input bin FASTA (remove contigs < 5kb)
2. FASTA index preparation
3. Align Hi-C reads to filtered FASTA OR use provided BAM
4. Scaffolding using YaHS
5. Cut scaffolds into fixed-length segments
6. Hi-C read realignment to scaffolded genome
7. MetaCC contact matrix generation
8. Contact matrix normalization
9. Scaffold mapping creation
10. Heatmap visualization
11. Quality metrics calculation

Dependencies:
- samtools
- YaHS
- MetaCC
- Matplotlib
- CheckM2
- pysam
"""

import argparse
import logging
import sys
import os
import subprocess
import warnings
import gzip
import pickle
import pandas as pd
import scipy.sparse as scisp
from pandas.errors import SettingWithCopyWarning
warnings.filterwarnings("ignore", category=SettingWithCopyWarning)

from raw_contact_both import ContactMatrix
from MetaCC.Script.normalized_contact import NormCCMap
from MetaCC.Script.utils import save_object, make_dir
from MetaCC.Script.normcc import normcc

class ApplicationException(Exception):
    def __init__(self, message):
        super(ApplicationException, self).__init__(message)

def filter_fasta_by_length(input_fasta, output_fasta, min_length=5000):
    """
    Filter FASTA file to remove contigs shorter than min_length
    
    Args:
        input_fasta (str): Path to input FASTA file
        output_fasta (str): Path to output filtered FASTA file
        min_length (int): Minimum contig length to keep
    
    Returns:
        int: Number of contigs retained
    """
    retained_count = 0
    current_seq = ""
    current_header = ""
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Process previous sequence
                if current_header and len(current_seq) >= min_length:
                    outfile.write(current_header + '\n')
                    outfile.write(current_seq + '\n')
                    retained_count += 1
                
                # Start new sequence
                current_header = line.strip()
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Don't forget the last sequence
        if current_header and len(current_seq) >= min_length:
            outfile.write(current_header + '\n')
            outfile.write(current_seq + '\n')
            retained_count += 1
    
    return retained_count

def cut_scaffolds_to_segments(scaffolded_fasta, output_fasta, segment_length):
    """
    Cut scaffolds into fixed-length segments for visualization
    
    Args:
        scaffolded_fasta (str): Path to scaffolded FASTA file
        output_fasta (str): Path to output segmented FASTA file
        segment_length (int): Length of each segment
    
    Returns:
        dict: Mapping of segment names to scaffold names
    """
    segment_to_scaffold = {}
    
    with open(scaffolded_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        current_scaffold = None
        current_seq = ""
        
        for line in infile:
            if line.startswith('>'):
                # Process previous scaffold
                if current_scaffold and current_seq:
                    # Cut into segments
                    for i in range(0, len(current_seq), segment_length):
                        segment = current_seq[i:i+segment_length]
                        if len(segment) > 0:  
                            segment_name = f"{current_scaffold}_seg_{i//segment_length + 1}"
                            outfile.write(f">{segment_name}\n{segment}\n")
                            segment_to_scaffold[segment_name] = current_scaffold
                
                # Start new scaffold
                current_scaffold = line[1:].strip().split()[0]
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Process last scaffold
        if current_scaffold and current_seq:
            for i in range(0, len(current_seq), segment_length):
                segment = current_seq[i:i+segment_length]
                if len(segment) > 0:  
                    segment_name = f"{current_scaffold}_seg_{i//segment_length + 1}"
                    outfile.write(f">{segment_name}\n{segment}\n")
                    segment_to_scaffold[segment_name] = current_scaffold
    
    return segment_to_scaffold

def create_scaffold_mapping(segment_to_scaffold_dict, output_path):
    """
    Create a mapping where each scaffold's segments are treated as contigs within that scaffold's 'bin'
    
    Args:
        segment_to_scaffold_dict (dict): Mapping of segment names to scaffold names
        output_path (str): Path to save the scaffold mapping
    
    Returns:
        dict: Mapping of segment names to scaffold names
    """
    # Save as both .p and .p.gz for compatibility
    save_object(output_path.replace('.gz', ''), segment_to_scaffold_dict)
    
    with gzip.open(output_path, 'wb') as f:
        pickle.dump(segment_to_scaffold_dict, f)
    
    return segment_to_scaffold_dict

def calculate_n50_l50(fasta_file):
    """
    Calculate N50 and L50 metrics for a FASTA file
    
    Args:
        fasta_file (str): Path to FASTA file
    
    Returns:
        tuple: (N50, L50)
    """
    lengths = []
    current_length = 0
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_length > 0:
                    lengths.append(current_length)
                current_length = 0
            else:
                current_length += len(line.strip())
    
    # Add the last sequence
    if current_length > 0:
        lengths.append(current_length)
    
    if not lengths:
        return 0, 0
    
    # Sort lengths in descending order
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    half_total = total_length / 2
    
    cumulative_length = 0
    for i, length in enumerate(lengths):
        cumulative_length += length
        if cumulative_length >= half_total:
            return length, i + 1  # N50, L50
    
    return 0, 0

def run_checkm2(fasta_file, output_dir, threads=1):
    """
    Run CheckM2 on a FASTA file
    
    Args:
        fasta_file (str): Path to FASTA file
        output_dir (str): Output directory for CheckM2
        threads (int): Number of threads
    
    Returns:
        tuple: (completeness, contamination)
    """
    checkm2_cmd = f"conda run -n checkm2 checkm2 predict --threads {threads} --input {fasta_file} --output-directory {output_dir}"
    
    try:
        result = subprocess.run(checkm2_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"CheckM2 failed: {result.stderr}")
            return None, None
        
        # Parse CheckM2 output
        quality_file = os.path.join(output_dir, "quality_report.tsv")
        if os.path.exists(quality_file):
            df = pd.read_csv(quality_file, sep='\t')
            if not df.empty:
                completeness = df['Completeness'].iloc[0]
                contamination = df['Contamination'].iloc[0]
                return completeness, contamination
    except Exception as e:
        logger.warning(f"CheckM2 error: {e}")
    
    return None, None

def ensure_bam_indexed(bam_file):
    """
    Create BAM index if it doesn't exist
    
    Args:
        bam_file (str): Path to BAM file
    
    Returns:
        bool: True if index exists or was created successfully
    """
    index_file = bam_file + ".bai"
    if not os.path.exists(index_file):
        logger.info(f"Creating BAM index for {bam_file}")
        exit_code = os.system(f"samtools index {bam_file}")
        if exit_code != 0:
            logger.error(f"Failed to index {bam_file}")
            return False
    return True

def create_coordinate_sorted_bam(input_bam, output_bam, threads=1):
    """
    Create coordinate-sorted BAM from input BAM
    
    Args:
        input_bam (str): Path to input BAM file
        output_bam (str): Path to output coordinate-sorted BAM file
        threads (int): Number of threads for sorting
    
    Returns:
        bool: True if sorting was successful
    """
    if os.path.exists(output_bam):
        logger.info(f"Coordinate-sorted BAM already exists: {output_bam}")
        return True
    
    logger.info(f"Creating coordinate-sorted BAM from {input_bam}")
    sort_cmd = f"samtools sort -@ {threads} {input_bam} -o {output_bam}"
    
    exit_code = os.system(sort_cmd)
    if exit_code != 0:
        logger.error(f"Failed to create coordinate-sorted BAM: {output_bam}")
        return False
    
    logger.info(f"Coordinate-sorted BAM created: {output_bam}")
    return True

def calculate_hic_enrichment_contigs(bam_file, filtered_fasta):
    """
    Calculate Hi-C enrichment ratio for original contigs (intra-contig / inter-contig contacts)
    
    Args:
        bam_file (str): Path to BAM file with Hi-C reads aligned to original contigs
        filtered_fasta (str): Path to original filtered FASTA file
    
    Returns:
        float: Hi-C enrichment ratio
    """
    try:
        import pysam
        
        # Ensure BAM is indexed
        if not ensure_bam_indexed(bam_file):
            logger.warning("Cannot index BAM file for contig enrichment calculation")
            return None
        
        # Read contig names from filtered FASTA
        original_contigs = set()
        with open(filtered_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_name = line[1:].strip().split()[0]
                    original_contigs.add(contig_name)
        
        logger.info(f"Found {len(original_contigs)} original contigs")
        
        # Analyze BAM file to count contacts
        within_contig_contacts = 0
        across_contig_contacts = 0
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch():
                # Remove proper pair requirement - just check if both reads are mapped and paired
                if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
                    read_contig = read.reference_name
                    mate_contig = read.next_reference_name
                    
                    # Both reads must map to contigs we're analyzing
                    if read_contig in original_contigs and mate_contig in original_contigs:
                        if read_contig == mate_contig:
                            within_contig_contacts += 1
                        else:
                            across_contig_contacts += 1
        
        # Calculate enrichment ratio
        if across_contig_contacts > 0:
            enrichment_ratio = within_contig_contacts / across_contig_contacts
        else:
            enrichment_ratio = float('inf') if within_contig_contacts > 0 else 0
        
        logger.info(f"Original contigs - Within: {within_contig_contacts}, Across: {across_contig_contacts}")
        return enrichment_ratio
    
    except Exception as e:
        logger.warning(f"Original contig Hi-C enrichment calculation error: {e}")
        return None


def calculate_hic_enrichment_scaffolds(bam_file, scaffolded_fasta):
    """
    Calculate Hi-C enrichment ratio for scaffolds (intra-scaffold / inter-scaffold contacts)
    
    Args:
        bam_file (str): Path to BAM file with Hi-C reads aligned to segmented scaffolds  
        scaffolded_fasta (str): Path to original scaffolded FASTA file
    
    Returns:
        float: Hi-C enrichment ratio
    """
    try:
        import pysam
        
        # Ensure BAM is indexed
        if not ensure_bam_indexed(bam_file):
            logger.warning("Cannot index BAM file for scaffold enrichment calculation")
            return None
        
        # Read scaffold names from scaffolded FASTA
        scaffolds = set()
        with open(scaffolded_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    scaffold_name = line[1:].strip().split()[0]
                    scaffolds.add(scaffold_name)
        
        logger.info(f"Found {len(scaffolds)} scaffolds")
        
        # Function to extract scaffold name from segment name
        def get_scaffold_from_segment(segment_name):
            if '_seg_' in segment_name:
                return segment_name.split('_seg_')[0]
            return segment_name
        
        # Analyze BAM file to count contacts
        within_scaffold_contacts = 0
        across_scaffold_contacts = 0
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch():
                # Remove proper pair requirement - just check if both reads are mapped and paired
                if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
                    # Get segment names from BAM
                    read_segment = read.reference_name
                    mate_segment = read.next_reference_name
                    
                    # Extract scaffold names from segment names
                    read_scaffold = get_scaffold_from_segment(read_segment)
                    mate_scaffold = get_scaffold_from_segment(mate_segment)
                    
                    # Both reads must map to scaffolds we're analyzing
                    if read_scaffold in scaffolds and mate_scaffold in scaffolds:
                        if read_scaffold == mate_scaffold:
                            within_scaffold_contacts += 1
                        else:
                            across_scaffold_contacts += 1
        
        # Calculate enrichment ratio
        if across_scaffold_contacts > 0:
            enrichment_ratio = within_scaffold_contacts / across_scaffold_contacts
        else:
            enrichment_ratio = float('inf') if within_scaffold_contacts > 0 else 0
        
        logger.info(f"Scaffolds - Within: {within_scaffold_contacts}, Across: {across_scaffold_contacts}")
        return enrichment_ratio
    
    except Exception as e:
        logger.warning(f"Scaffold Hi-C enrichment calculation error: {e}")
        return None


def calculating_metrics(filtered_fasta, genome_path, contact_matrix, segment_to_scaffold, args, output_dir):
    """
    Calculate scaffolding quality metrics using direct BAM analysis
    """
    logger.info("Calculating scaffolding quality metrics...")
    
    # Calculate N50 and L50 for original and scaffolded genomes
    original_n50, original_l50 = calculate_n50_l50(filtered_fasta)
    scaffolded_n50, scaffolded_l50 = calculate_n50_l50(genome_path)
    
    # Count total contigs and scaffolds
    original_contig_count = 0
    with open(filtered_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                original_contig_count += 1
    
    scaffolded_scaffold_count = 0
    with open(genome_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                scaffolded_scaffold_count += 1
    
    # Run CheckM2 on both original and scaffolded genomes
    checkm2_original_dir = os.path.join(output_dir, "checkm2_original")
    checkm2_scaffolded_dir = os.path.join(output_dir, "checkm2_scaffolded")
    make_dir(checkm2_original_dir)
    make_dir(checkm2_scaffolded_dir)
    
    threads = int(args.t) if args.t else 1
    
    logger.info("Running CheckM2 on original bin...")
    original_completeness, original_contamination = run_checkm2(filtered_fasta, checkm2_original_dir, threads)
    
    logger.info("Running CheckM2 on scaffolded genome...")
    scaffolded_completeness, scaffolded_contamination = run_checkm2(genome_path, checkm2_scaffolded_dir, threads)
    
    # Calculate Hi-C enrichment ratios using direct BAM analysis
    logger.info("Calculating Hi-C enrichment ratios...")
    
    # Determine which BAM to use for original contig enrichment
    original_hic_enrichment = None
    scaffolded_hic_enrichment = None
    
    # For original contigs: prioritize coordinate-sorted versions
    coord_sorted_original = os.path.join(output_dir, "sorted_for_enrichment.bam")
    initial_coord_sorted = os.path.join(output_dir, "initial_sorted_for_enrichment.bam")
    initial_alignment_bam = os.path.join(output_dir, "initial_alignment", "sorted_map.bam")
    yahs_sorted_bam = os.path.join(output_dir, "sorted_for_yahs.bam")
    
    if os.path.exists(coord_sorted_original):
        logger.info("Using coordinate-sorted BAM (from provided BAM) for original contig enrichment")
        original_hic_enrichment = calculate_hic_enrichment_contigs(coord_sorted_original, filtered_fasta)
    elif os.path.exists(initial_coord_sorted):
        logger.info("Using coordinate-sorted initial alignment BAM for original contig enrichment")
        original_hic_enrichment = calculate_hic_enrichment_contigs(initial_coord_sorted, filtered_fasta)
    elif os.path.exists(initial_alignment_bam):
        logger.info("Creating coordinate-sorted BAM from initial alignment for enrichment analysis")
        if create_coordinate_sorted_bam(initial_alignment_bam, initial_coord_sorted, threads):
            original_hic_enrichment = calculate_hic_enrichment_contigs(initial_coord_sorted, filtered_fasta)
        else:
            logger.warning("Failed to create coordinate-sorted BAM from initial alignment")
    elif os.path.exists(yahs_sorted_bam):
        logger.info("Creating coordinate-sorted BAM from YaHS BAM for enrichment analysis")
        if create_coordinate_sorted_bam(yahs_sorted_bam, coord_sorted_original, threads):
            original_hic_enrichment = calculate_hic_enrichment_contigs(coord_sorted_original, filtered_fasta)
        else:
            logger.warning("Failed to create coordinate-sorted BAM from YaHS BAM")
    else:
        logger.warning("No suitable BAM found for original contig enrichment")
    
    # For scaffolds: prioritize coordinate-sorted realigned BAM
    realigned_coord_bam = os.path.join(output_dir, "realign_sorted_for_enrichment.bam")
    realigned_bam = os.path.join(output_dir, "realign_to_scaffold", "sorted_map.bam")
    
    if os.path.exists(realigned_coord_bam):
        logger.info("Using coordinate-sorted realigned BAM for scaffold enrichment")
        scaffolded_hic_enrichment = calculate_hic_enrichment_scaffolds(realigned_coord_bam, genome_path)
    elif os.path.exists(realigned_bam):
        logger.info("Creating coordinate-sorted BAM from realigned BAM for scaffold enrichment")
        if create_coordinate_sorted_bam(realigned_bam, realigned_coord_bam, threads):
            scaffolded_hic_enrichment = calculate_hic_enrichment_scaffolds(realigned_coord_bam, genome_path)
        else:
            logger.warning("Failed to create coordinate-sorted BAM from realigned BAM")
    else:
        logger.warning("No suitable BAM found for scaffold enrichment")
    
    # Compile metrics
    metrics = {
        "Original_N50": original_n50,
        "Original_L50": original_l50,
        "Original_contig_count": original_contig_count,
        "Scaffolded_N50": scaffolded_n50,
        "Scaffolded_L50": scaffolded_l50,
        "Scaffolded_scaffold_count": scaffolded_scaffold_count,
        "Original_completeness": original_completeness,
        "Original_contamination": original_contamination,
        "Scaffolded_completeness": scaffolded_completeness,
        "Scaffolded_contamination": scaffolded_contamination,
        "Original_HiC_enrichment_ratio": original_hic_enrichment,
        "Scaffolded_HiC_enrichment_ratio": scaffolded_hic_enrichment
    }
    
    return metrics

def save_metrics(metrics, output_file):
    """
    Save metrics to a text file
    
    Args:
        metrics (dict): Dictionary of metrics
        output_file (str): Path to output file
    """
    with open(output_file, 'w') as f:
        f.write("Scaffolding Quality Metrics\n")
        f.write("==========================\n\n")
        
        for key, value in metrics.items():
            if value is not None:
                if isinstance(value, float):
                    f.write(f"{key}: {value:.4f}\n")
                else:
                    f.write(f"{key}: {value}\n")
            else:
                f.write(f"{key}: N/A\n")

if __name__ == '__main__':

    def make_dir(path, exist_ok=True):
        if os.path.isfile(path):
            raise IOError('Output path already exists and is a file!')
        if os.path.isdir(path):
            if not exist_ok:
                raise IOError('Output directory already exists!')
        else:
            os.makedirs(path)

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', help='Base metahit path')
    parser.add_argument('--fasta', help='Reference fasta sequence of bin for scaffolding')
    parser.add_argument('--bam', help='Input bam file of Hi-C reads mapped to shotgun assembly (optional if using --hic1 and --hic2)')
    parser.add_argument('--enzyme', help='List of restriction enzymes, separated by comma')
    parser.add_argument('--outdir', help='Output directory path')
    parser.add_argument('--hic1', help='Hi-C forward fastq')
    parser.add_argument('--hic2', help='Hi-C reverse fastq')
    parser.add_argument('-t', help='Threads')
    parser.add_argument('-m', help='Memory')
    parser.add_argument('-r', help='Resolution (default 10kb)')

    args = parser.parse_args()

    try:
        make_dir(args.outdir)
    except IOError:
        print('Error: cannot find output directory or it already exists')
        sys.exit(1)

    logging.captureWarnings(True)
    logger = logging.getLogger('main')
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    log_path = os.path.join(args.outdir, 'scaffolding.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:
        # Step 1: Filter input FASTA to remove short contigs
        filtered_fasta = os.path.join(args.outdir, "filtered_bin.fa")
        MIN_CONTIG_LEN = 5000
        
        logger.info(f"Filtering input FASTA to remove contigs < {MIN_CONTIG_LEN} bp...")
        retained_contigs = filter_fasta_by_length(args.fasta, filtered_fasta, MIN_CONTIG_LEN)
        logger.info(f"Retained {retained_contigs} contigs after filtering")
        
        if retained_contigs == 0:
            logger.error("No contigs retained after filtering. Exiting.")
            sys.exit(1)

        # Step 2: Ensure .fai index exists for filtered FASTA
        fai_path = filtered_fasta + ".fai"
        if not os.path.exists(fai_path):
            logger.info(f"Creating FASTA index for filtered bin...")
            exit_code = os.system(f"samtools faidx {filtered_fasta}")
            if exit_code != 0:
                logger.error("samtools faidx failed. Cannot proceed.")
                sys.exit(1)
            logger.info("FASTA index created.")
        else:
            logger.info("FASTA index already exists.")

        # Step 3: Align Hi-C reads to filtered FASTA or use provided BAM
        if args.bam:
            # Use provided BAM file
            logger.info(f"Using provided BAM file: {args.bam}")
            # Create name-sorted BAM for YaHS
            sorted_bam = os.path.join(args.outdir, "sorted_for_yahs.bam")
            sort_cmd = f"samtools sort -n {args.bam} -o {sorted_bam}"
            logger.info("Creating name-sorted BAM for YaHS: " + sort_cmd)
            
            if os.system(sort_cmd) != 0:
                logger.error("BAM sorting failed.")
                sys.exit(1)
        else:
            # Align Hi-C reads to filtered FASTA
            logger.info("No BAM provided. Aligning Hi-C reads to filtered FASTA...")
            if not args.hic1 or not args.hic2:
                logger.error("Hi-C reads (--hic1 and --hic2) are required when no BAM is provided.")
                sys.exit(1)
            
            # Create alignment directory
            initial_alignment_dir = os.path.join(args.outdir, "initial_alignment")
            make_dir(initial_alignment_dir)
            
            # Run alignment to filtered FASTA
            alignment_cmd = (
                f"bash {os.path.dirname(os.path.abspath(__file__))}/alignment.sh "
                f"-p {args.p} "
                f"-r {filtered_fasta} "
                f"-1 {args.hic1} "
                f"-2 {args.hic2} "
                f"-o {initial_alignment_dir} "
                f"-t {args.t}"
            )
            
            logger.info("Running initial alignment to filtered FASTA:\n" + alignment_cmd)
            if os.system(alignment_cmd) != 0:
                logger.error("Initial alignment failed.")
                sys.exit(1)
            
            # Use the sorted BAM from alignment
            sorted_bam = os.path.join(initial_alignment_dir, "sorted_map.bam")
            
            # Verify alignment output
            if not os.path.exists(sorted_bam):
                logger.error(f"Aligned BAM not found at {sorted_bam}")
                sys.exit(1)
                
            logger.info("Initial alignment completed successfully.")

        # Step 4: Run YaHS on filtered FASTA
        yahs_prefix = os.path.join(args.outdir, "yahs", "scaffold")
        make_dir(os.path.dirname(yahs_prefix))
        
        yahs_cmd = f"yahs {filtered_fasta} {sorted_bam} -o {yahs_prefix}"
        logger.info("Running YaHS: " + yahs_cmd)
        
        if os.system(yahs_cmd) != 0:
            logger.error("YaHS failed.")
            sys.exit(1)
        
        logger.info("YaHS completed successfully.")

        genome_path = yahs_prefix + "_scaffolds_final.fa"
        
        # Verify YaHS output exists
        if not os.path.exists(genome_path):
            logger.error(f"YaHS output not found at {genome_path}")
            sys.exit(1)

        # Step 5: Cut scaffolds into fixed-length segments
        segment_length = int(args.r) if args.r else 10000  
        logger.info(f"Using segment length: {segment_length} bp")
        
        segmented_fasta = os.path.join(args.outdir, "segmented_scaffolds.fa")
        segment_to_scaffold = cut_scaffolds_to_segments(genome_path, segmented_fasta, segment_length)
        
        logger.info(f"Created {len(segment_to_scaffold)} segments from scaffolds")

        # Step 6: Realign Hi-C reads to scaffolded genome
        realigned_dir = os.path.join(args.outdir, "realign_to_scaffold")
        make_dir(realigned_dir)
        
        # Use the alignment script in the same directory as scaffolding.py
        alignment_cmd = (
            f"bash {os.path.dirname(os.path.abspath(__file__))}/alignment.sh "
            f"-p {args.p} "
            f"-r {segmented_fasta} "  # Use segmented FASTA for alignment
            f"-1 {args.hic1} "
            f"-2 {args.hic2} "
            f"-o {realigned_dir} "
            f"-t {args.t}"
        )
        sorted_bam_path = os.path.join(realigned_dir, "sorted_map.bam")

        logger.info("Running realignment to segmented scaffolds:\n" + alignment_cmd)
        if os.system(alignment_cmd) != 0:
            logger.error("Realignment to segmented scaffolds failed.")
            sys.exit(1)
        
        # Verify realignment output
        if not os.path.exists(sorted_bam_path):
            logger.error(f"Realigned BAM not found at {sorted_bam_path}")
            sys.exit(1)
            
        logger.info("Realignment to segmented scaffolds completed successfully.")
            
        # Step 7: Build MetaCC ContactMatrix using segmented scaffolds
        metacc_folder = os.path.join(args.outdir, "metacc")
        make_dir(metacc_folder)
        
        # Create tmp folder structure that ContactMatrix expects
        metacc_temp_folder = os.path.join(metacc_folder, 'tmp')
        make_dir(metacc_temp_folder)
    
        logger.info("Generating MetaCC contact map for segmented scaffolds...")
    
        enzymes = args.enzyme.split(",") if args.enzyme else ["Sau3AI", "MluCI"]
        
        contact_matrix = ContactMatrix(
            sorted_bam_path,  # Use realigned BAM
            enzymes,
            segmented_fasta,  # Use segmented scaffolds
            args.outdir,
            metacc_folder,
            None,             # bin3c_folder 
            min_insert=0,
            bin_size=None,
            min_mapq_metacc=30,
            min_len_metacc=1000,
            min_match_metacc=30,
            min_signal_metacc=2,
            min_mapq_bin3c=60,
            min_len_bin3c=1000,
            min_match_bin3c=10,
            min_signal_bin3c=5
        )
    
        # Save contact matrix with both extensions for compatibility
        cm_path = os.path.join(metacc_folder, "contact_map.p")
        cm_path_gz = os.path.join(metacc_folder, "contact_map.p.gz")
        save_object(cm_path, contact_matrix)
        
        with gzip.open(cm_path_gz, 'wb') as f:
            pickle.dump(contact_matrix, f)
            
        logger.info(f"MetaCC contact map saved to: {cm_path}")
        
        # Step 8: Normalize contact matrix with MetaCC
        logger.info('Begin normalizing raw contacts by NormCC for scaffolded matrix')
        
        contig_file = os.path.join(metacc_temp_folder, 'contig_info_metacc.csv')
        
        # Verify contig file exists (should be created by ContactMatrix in tmp folder)
        if not os.path.exists(contig_file):
            logger.error(f"Contig info file not found at {contig_file}")
            # Check what files are actually in the tmp folder
            logger.error(f"Files in tmp folder: {os.listdir(metacc_temp_folder) if os.path.exists(metacc_temp_folder) else 'tmp folder does not exist'}")
            # Check what files are in metacc folder
            logger.error(f"Files in metacc folder: {os.listdir(metacc_folder) if os.path.exists(metacc_folder) else 'metacc folder does not exist'}")
            sys.exit(1)
            
        norm_result = normcc(contig_file)
        contig_info_df = pd.read_csv(contig_file)
        
        hzmap = NormCCMap(
            metacc_folder,
            contig_info_df,
            contact_matrix.seq_map_metacc,
            norm_result,
            thres=0.05
        )
        
        logger.info('NormCC normalization finished')
        
        # Save normalized matrix
        normalized_matrix_path = os.path.join(metacc_folder, 'Normalized_contact_matrix.npz')
        scisp.save_npz(normalized_matrix_path, hzmap.seq_map.tocsr())
        
        normalized_obj_path = os.path.join(metacc_folder, 'NormCC_normalized_contact.gz')
        with gzip.open(normalized_obj_path, 'wb') as f:
            pickle.dump(hzmap, f)
        
        logger.info('MetaCC Normalization results have been saved')

        # Step 9: Create scaffold mapping for visualization
        logger.info('Creating scaffold mapping for heatmap visualization...')
        final_bins_path = os.path.join(args.outdir, "final_bins.p.gz")
        scaffold_mapping = create_scaffold_mapping(segment_to_scaffold, final_bins_path)
        logger.info(f'Scaffold mapping created with {len(segment_to_scaffold)} segments across {len(set(segment_to_scaffold.values()))} scaffolds')

        # Step 10: Generate heatmap
        heatmap_outdir = os.path.join(args.outdir, "figures")
        make_dir(heatmap_outdir)
        
        # Verify all required files exist
        required_files = [
            normalized_matrix_path,
            cm_path_gz,
            final_bins_path
        ]
        
        for file_path in required_files:
            if not os.path.exists(file_path):
                logger.error(f"Required file not found: {file_path}")
                sys.exit(1)
        
        heatmap_cmd = (
            f"python {args.p}/8_scaffolding/scripts/heatmap.py "
            f"--contact-map {normalized_matrix_path} "
            f"--ORDER {cm_path_gz} "
            f"--BIN {final_bins_path} "
            f"--OUTDIR {heatmap_outdir}"
        )
        
        logger.info("Running heatmap generation:\n" + heatmap_cmd)
        if os.system(heatmap_cmd) != 0:
            logger.error("Heatmap generation failed.")
            sys.exit(1)

        # Step 11: Calculate quality metrics
        metrics = calculating_metrics(filtered_fasta, genome_path, contact_matrix, segment_to_scaffold, args, args.outdir)
        
        # Save metrics
        metrics_file = os.path.join(args.outdir, "scaffolding_metrics.txt")
        save_metrics(metrics, metrics_file)
        logger.info(f"Scaffolding metrics saved to: {metrics_file}")
        
        logger.info("Scaffolding pipeline completed successfully!")
        logger.info(f"Scaffolded genome: {genome_path}")
        logger.info(f"Segmented scaffolds: {segmented_fasta}")
        logger.info(f"Contact matrix: {normalized_matrix_path}")
        logger.info(f"Heatmap: {os.path.join(heatmap_outdir, 'heatmap.png')}")
        logger.info(f"Metrics: {metrics_file}")

    except ApplicationException as e:
        logger.error(f"ApplicationException Error: {str(e)}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        sys.exit(1)