#!/usr/bin/env python
import sys
import textwrap

records = {}
name = None
seq = []

with open(sys.argv[1]) as handle:
    for line in handle:
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seq)
            name = line.strip()
            seq = []
        else:
            seq.append(line.strip())

if name is not None:
    records[name] = "".join(seq)

for header in sorted(records, key=lambda key: len(records[key]), reverse=True):
    print(header)
    print(textwrap.fill(records[header], 100, break_on_hyphens=False))
