#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import os
import subprocess
import numpy as np
from scipy.stats import norm


def parse_args():
    # Get the script's directory
    script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    # Go two directories up to get to the metahit path
    default_metahit_path = os.path.abspath(os.path.join(script_dir, '../..'))
    
    parser = argparse.ArgumentParser(
        description="EM-based Hi-C/WGS split and reassembly"
    )
    parser.add_argument('--bin', required=True, help='Binning result directory')
    parser.add_argument('--hic1', required=True, help='Hi-C library forward FASTQ')
    parser.add_argument('--hic2', required=True, help='Hi-C library reverse FASTQ')
    parser.add_argument('--sg1', required=True, help='Shotgun forward FASTQ')
    parser.add_argument('--sg2', required=True, help='Shotgun reverse FASTQ')
    parser.add_argument('--bam', required=True, help='Hi-C BAM mapped to assembly')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('-p', '--metahit_path', default=default_metahit_path,
                        help='Path to metahit folder (default: two directories up from script)')
    parser.add_argument('-t', '--threads', default='4', help='Number of threads')
    parser.add_argument('-m', '--memory', default='24', help='Memory (GB)')
    parser.add_argument('-k', '--top_k', type=int, default=100,
                        help='Number of top contigs to extract insert sizes from')
    parser.add_argument('-c', '--contig_sort_list', required=True,
                        help='File listing contig names (one per line) sorted by length (desc)')
    return parser.parse_args()

def run_command(cmd, check=True):
    """Run shell command with error handling"""
    print(f"[CMD] {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=check, capture_output=True, text=True)
        if result.stdout:
            print(f"[STDOUT] {result.stdout.strip()}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}")
        print(f"[STDERR] {e.stderr}")
        raise

def sample_insert_sizes(bam_path, output_file, top_k, contig_sort_file):
    """
    Extract insert sizes for all reads mapped to the top_k contigs specified in contig_sort_file.
    contig_sort_file: text file with one contig name per line, sorted by length descending.
    """
    print(f"[INFO] Sampling insert sizes from top {top_k} contigs...")
    
    contig_list = os.path.join(os.path.dirname(output_file), 'top_contigs.txt')
    run_command(f"head -n {top_k} {contig_sort_file} > {contig_list}")
    
    # More robust approach for extracting insert sizes
    extract_cmd = (
        f"samtools view {bam_path} | "
        f"grep -F -f {contig_list} | "
        f"awk '$9 != 0 && $9 > -5000 && $9 < 5000 {{if($9 < 0) print -$9; else print $9}}' > {output_file}"
    )
    run_command(extract_cmd)
    
    # Check if file is empty
    if os.path.getsize(output_file) == 0:
        print("[WARNING] No insert sizes found. Using default values.")
        with open(output_file, 'w') as f:
            f.write("300\n400\n500\n")

def init_params(data, init_frac=0.8):
    """
    Initialize EM parameters by splitting data at the init_frac quantile.
    init_frac=0.8 means lower 80% as non-chimeric.
    """
    if len(data) == 0:
        print("[WARNING] No insert size data available. Using default parameters.")
        return 400, 2000, 100, 500, 0.8, 0.2
    
    q = np.percentile(data, init_frac * 100)
    lower = data[data <= q]   # ~80% non-chimeric insert sizes
    upper = data[data > q]    # ~20% chimeric insert sizes
    
    if len(lower) == 0:
        lower = data[:max(1, len(data)//2)]
    if len(upper) == 0:
        upper = data[max(1, len(data)//2):]
    
    mu_N = lower.mean()
    mu_C = upper.mean()
    sigma_N = max(lower.std(ddof=1), 1.0)
    sigma_C = max(upper.std(ddof=1), 1.0)
    pi_N = len(lower) / len(data)
    pi_C = 1.0 - pi_N
    return mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C

def em_mix(data, mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C, tol=1e-2, max_iter=100):
    """Expectation-Maximization for mixture model of insert sizes"""
    if len(data) == 0:
        print("[WARNING] No data for EM algorithm. Returning initial parameters.")
        return mu_C, mu_N, sigma_C, sigma_N, pi_C, pi_N
    
    logL = None
    for iteration in range(max_iter):
        # E-step: Compute responsibilities
        wC = pi_C * norm.pdf(data, mu_C, sigma_C)
        wN = pi_N * norm.pdf(data, mu_N, sigma_N)
        total = wC + wN
        
        # Avoid division by zero
        total = np.maximum(total, 1e-10)
        
        gamma_C = wC / total
        gamma_N = wN / total
        
        # M-step: Update parameters
        sum_gamma_C = np.sum(gamma_C)
        sum_gamma_N = np.sum(gamma_N)
        
        if sum_gamma_C > 0:
            mu_C = np.sum(gamma_C * data) / sum_gamma_C
            sigma_C = max(np.sqrt(np.sum(gamma_C * (data - mu_C)**2) / sum_gamma_C), 1e-6)
        
        if sum_gamma_N > 0:
            mu_N = np.sum(gamma_N * data) / sum_gamma_N
            sigma_N = max(np.sqrt(np.sum(gamma_N * (data - mu_N)**2) / sum_gamma_N), 1e-6)
        
        pi_C = gamma_C.mean()
        pi_N = gamma_N.mean()
        
        new_logL = np.sum(np.log(pi_C * norm.pdf(data, mu_C, sigma_C) +
                              pi_N * norm.pdf(data, mu_N, sigma_N)))
        
        if logL is not None and abs(new_logL - logL) < tol:
            print(f"[INFO] EM converged after {iteration+1} iterations")
            break
        logL = new_logL
        print(f'Iter {iteration+1}: logL={new_logL:.2f}, pi_C={pi_C:.2f}, pi_N={pi_N:.2f}, mu_C={mu_C:.1f}, mu_N={mu_N:.1f}')
    
    return mu_C, mu_N, sigma_C, sigma_N, pi_C, pi_N

def count_intra_inter_pairs(bam_path):
    """
    Count read pairs mapped to the same contig (n) and different contigs (m).
    Only counts primary alignments to avoid double-counting.
    """
    print("[INFO] Counting intra/inter-contig pairs...")
    
    n = 0  # intra-contig
    m = 0  # inter-contig
    
    # Filter out secondary (0x100) and supplementary (0x800) alignments
    # -f 1: only paired reads
    # -F 0x904: exclude unmapped (0x4), secondary (0x100), supplementary (0x800)
    cmd = f"samtools view -f 1 -F 0x904 {bam_path}"
    result = run_command(cmd, check=False)
    
    # Use a set to track unique read names to ensure each pair is counted only once
    seen = set()
    
    for line in result.stdout.split('\n'):
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 7:
            continue
            
        qname = fields[0]  # read name
        rname = fields[2]  # reference name (contig)
        mate_rname = fields[6]  # mate reference name
        
        # Skip if we've already processed this read pair
        if qname in seen:
            continue
        seen.add(qname)
        
        # Count intra vs inter-contig pairs
        if mate_rname == '=' or mate_rname == rname:
            n += 1
        else:
            m += 1
    
    print(f"[INFO] Intra-contig pairs: {n}, Inter-contig pairs: {m}")
    return n, m

def main():
    args = parse_args()
    
    # Check that inputs are .gz files
    for name, path in {
        "--hic1": args.hic1,
        "--hic2": args.hic2,
        "--sg1": args.sg1,
        "--sg2": args.sg2
    }.items():
        if not path.endswith(".gz"):
            raise ValueError(f"Error: {name} must be a gzip-compressed (.fastq.gz) file. Got: {path}")

    print(f"[INFO] Creating output directory: {args.outdir}")
    os.makedirs(args.outdir, exist_ok=True)

    # Insert size analysis
    print("[INFO] Starting insert size analysis...")
    insert_file = os.path.join(args.outdir, 'insert_size.txt')
    sample_insert_sizes(args.bam, insert_file, args.top_k, args.contig_sort_list)

    values = []
    with open(insert_file) as f:
        for line in f:
            line = line.strip()
            if line and line.replace('.', '').isdigit():
                val = float(line)
                if val > 0 and val < 5000:
                    values.append(val)
    
    data = np.array(values)
    print(f"[INFO] Found {len(data)} valid insert size values")

    mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C = init_params(data)
    mu_C, mu_N, sigma_C, sigma_N, pi_C, pi_N = em_mix(
        data, mu_N, mu_C, sigma_N, sigma_C, pi_N, pi_C
    )
    
    # Save estimated pi_C to file
    pi_file = os.path.join(args.outdir, "mixing_proportion.txt")
    with open(pi_file, "w") as f:
        f.write(f"{pi_C}\n")
    print(f"[INFO] Estimated pi_C saved to {pi_file}")
    
    # Calculate and save long-range chimeric read ratio
    n, m = count_intra_inter_pairs(args.bam)
    long_range_ratio = (n * pi_C + m) / (n + m) if (n + m) > 0 else float('nan')
    longrange_file = os.path.join(args.outdir, "long_range_ratio.txt")
    with open(longrange_file, "w") as f:
        f.write(f"{long_range_ratio}\n")
    print(f"[INFO] Long-range chimeric read ratio saved to {longrange_file}")
    
    ratio_3d = m / n if n > 0 else float('inf')
    ratio3d_file = os.path.join(args.outdir, "3d_ratio.txt")
    with open(ratio3d_file, "w") as f:
        f.write(f"{ratio_3d:.4f}\n")
    print(f"[INFO] 3D ratio saved to {ratio3d_file}")

    # Calculate cutoff more robustly
    try:
        th_C = mu_C + norm.ppf(0.05) * sigma_C
        th_N = mu_N + norm.ppf(0.95) * sigma_N
        
        if th_N <= 0:
            cutoff = max(th_C, 100)  # Minimum cutoff of 100
        elif th_C <= 0:
            cutoff = max(th_N, 100)
        else:
            cutoff = max(min(th_C, th_N), 100)
    except:
        cutoff = 500  # Default cutoff
        
    print(f"[INFO] Insert size cutoff: {cutoff}")
    
    cutoff_file = os.path.join(args.outdir, "insert_size_cutoff.txt")
    with open(cutoff_file, "w") as f:
        f.write(f"{cutoff:.2f}\n")
    print(f"[INFO] Insert size cutoff saved to {cutoff_file}")
    
    # Extract shotgun reads from Hi-C
    print("[INFO] Extracting shotgun reads from Hi-C data...")
    readnames = os.path.join(args.outdir, 'readname_sg_in_hic.txt')
    extract_sg_cmd = (
        f"samtools view {args.bam} | "
        f"awk '($9 >= -{cutoff} && $9 <= {cutoff} && $9 != 0)' | "
        f"cut -f1 | sort -u > {readnames}"
    )
    run_command(extract_sg_cmd)

    sg_in_hic_fwd = os.path.join(args.outdir, 'sg_in_hic.forward.fastq.gz')
    sg_in_hic_rev = os.path.join(args.outdir, 'sg_in_hic.reverse.fastq.gz')
    combined_sg_fwd = os.path.join(args.outdir, 'new_sg_forward.fastq.gz')
    combined_sg_rev = os.path.join(args.outdir, 'new_sg_reverse.fastq.gz')

    run_command(f"seqtk subseq {args.hic1} {readnames} | gzip > {sg_in_hic_fwd}")
    run_command(f"seqtk subseq {args.hic2} {readnames} | gzip > {sg_in_hic_rev}")
    run_command(f"cat {sg_in_hic_fwd} {args.sg1} > {combined_sg_fwd}")
    run_command(f"cat {sg_in_hic_rev} {args.sg2} > {combined_sg_rev}")

    # Save non-shotgun-like Hi-C readnames
    print("[INFO] Extracting non-shotgun-like Hi-C reads...")
    all_readnames = os.path.join(args.outdir, 'all_hic_readnames.txt')
    non_shotgun_readnames = os.path.join(args.outdir, 'readname_non_sg_in_hic.txt')
    
    # Extract all read names from BAM
    run_command(f"samtools view {args.bam} | cut -f1 | sort -u > {all_readnames}")
    
    # Get non-sg-like read names using set difference
    run_command(f"comm -23 {all_readnames} {readnames} > {non_shotgun_readnames}")
    
    # Extract non-sg-like Hi-C reads
    non_sg_fwd = os.path.join(args.outdir, 'non_sg_in_hic.forward.fastq.gz')
    non_sg_rev = os.path.join(args.outdir, 'non_sg_in_hic.reverse.fastq.gz')
    run_command(f"seqtk subseq {args.hic1} {non_shotgun_readnames} | gzip > {non_sg_fwd}")
    run_command(f"seqtk subseq {args.hic2} {non_shotgun_readnames} | gzip > {non_sg_rev}")

    # Run reassembly (which now includes unmapped read extraction and assembly)
    print("[INFO] Starting bin reassembly and unmapped read processing...")
    reassemble_cmd = (
        f"bash {args.metahit_path}/7_reassembly/scripts/reassemble_bins.sh "
        f"-b {args.bin} -o {args.outdir} -1 {combined_sg_fwd} -2 {combined_sg_rev} "
        f"-m {args.memory} -t {args.threads} --parallel"
    )
    run_command(reassemble_cmd)
    
    print("[INFO] Reassembly pipeline completed successfully!")
    print(f"[INFO] Results:")
    print(f"  - Reassembled bins: {args.outdir}/reassembled_bins/")
    print(f"  - Unmapped contigs: {args.outdir}/unmapped_assembly/unmapped_contigs.fasta")
    print(f"  - Unmapped reads: {args.outdir}/unmapped_shotgun_1.fastq, {args.outdir}/unmapped_shotgun_2.fastq")

if __name__ == '__main__':
    main()