#!/usr/bin/env python
from __future__ import print_function
import sys

min_length = int(sys.argv[1])
current_contig = ""
current_sequence = ""

for line in open(sys.argv[2]):
    if line.startswith(">"):
        # Process the previous contig if there was one
        if current_contig and current_sequence:
            if len(current_sequence) >= min_length:
                print(current_contig)
                print(current_sequence)
        
        # Start a new contig
        current_contig = line.strip()
        current_sequence = ""
    else:
        current_sequence += line.strip()

# Process the last contig
if current_contig and current_sequence:
    if len(current_sequence) >= min_length:
        print(current_contig)
        print(current_sequence)