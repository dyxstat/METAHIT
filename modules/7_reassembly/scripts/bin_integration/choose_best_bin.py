#!/usr/bin/env python3
from __future__ import print_function
import sys
import pandas as pd

# This script takes in the reassembled_bins.stats file of the binning module and chooses the best possible
# scenario for each bin. Each bin should have 3 versions of itself: the original bin (.orig), the
# "strict" reassembled bin (.strict), and the "permissive" reassembled bin (.permissive). This program
# decides which one of the three is better for each bin.

if len(sys.argv) < 4:
    print("Usage: choose_best_bin.py <reassembled_bins.stats> <min_completion> <max_contamination>")
    sys.exit(1)

# Read stats file with pandas
try:
    stats_df = pd.read_csv(sys.argv[1], sep='\t')
    
    # Set minimum completion and maximum contamination
    min_completion = float(sys.argv[2])
    max_contamination = float(sys.argv[3])
    
    # Dictionary to store the best version of each bin
    best_bins = {}
    
    # Process each bin
    for _, row in stats_df.iterrows():
        # Skip header
        if "completeness" in str(row[0]).lower():
            continue
            
        bin_name = row.iloc[0]
        completeness = float(row.iloc[1])
        contamination = float(row.iloc[2])
        n50 = int(row.iloc[5])
        
        # Don't consider bins below min completion or above max contamination
        if completeness < min_completion or contamination > max_contamination:
            continue
            
        # Extract base bin name and style
        base_bin_name = ".".join(bin_name.split(".")[:-1])
        style = bin_name.split(".")[-1]
        
        # Calculate score (completeness + (100-contamination)*5)
        score = completeness + 5 * (100 - contamination)
        
        # Update best bin if this version is better
        if base_bin_name not in best_bins:
            best_bins[base_bin_name] = [style, score, n50]
        else:
            current_score = best_bins[base_bin_name][1]
            current_n50 = best_bins[base_bin_name][2]
            
            # Update if score is better, or if equal score but better N50
            if score > current_score or (score == current_score and n50 > current_n50):
                best_bins[base_bin_name] = [style, score, n50]
    
    # Output the best version of each bin
    for bin_name in best_bins:
        print(f"{bin_name}.{best_bins[bin_name][0]}")
        
except Exception as e:
    sys.stderr.write(f"Error processing stats file: {str(e)}\n")
    sys.exit(1)