#!/usr/bin/env python3
from __future__ import print_function
import sys
import pandas as pd

# This script summarizes the statistics of each bin by parsing 
# the quality_report.tsv file of the CheckM2 output

if len(sys.argv) == 3: 
    binner = sys.argv[2]
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner")
elif len(sys.argv) == 4:
    source = {}
    for line in open(sys.argv[3]):
        cut = line.strip().split("\t")
        source[cut[0]] = cut[7]
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner")
else:
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize")

try:
    # Read the quality report as a pandas DataFrame
    quality_df = pd.read_csv(sys.argv[1], sep='\t')
    
    # Process each bin
    for _, row in quality_df.iterrows():
        name = row['Name']
        completeness = row['Completeness']
        contamination = row['Contamination']
        gc_content = row['GC_Content']
        n50 = row['Contig_N50']
        size = row['Genome_Size']
        lineage = "Unknown"  # No Taxonomy available, so just use "Unknown"

        # Format output
        if len(sys.argv) == 3:
            print("\t".join([name, str(completeness)[:5],
                             str(contamination)[:5], str(gc_content)[:5],
                             lineage, str(n50), str(size), binner]))
        elif len(sys.argv) == 4:
            print("\t".join([name, str(completeness)[:5],
                             str(contamination)[:5], str(gc_content)[:5],
                             lineage, str(n50), str(size), source.get(name, "Unknown")]))
        else:
            print("\t".join([name, str(completeness)[:5],
                             str(contamination)[:5], str(gc_content)[:5],
                             lineage, str(n50), str(size)]))
except Exception as e:
    sys.stderr.write(f"Error processing CheckM2 output: {str(e)}\n")
    sys.exit(1)
