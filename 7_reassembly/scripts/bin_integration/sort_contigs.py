#!/usr/bin/env python
from __future__ import print_function
import sys
import textwrap

dic = {}
tmp_contig = ""
name = ""  

for line in open(sys.argv[1]):
    if line[0] == ">":
        if tmp_contig != "" and name != "":  
            dic[name] = tmp_contig
            tmp_contig = ""
        name = line.strip()
    else:
        tmp_contig += line.strip()

# Add the last contig
if tmp_contig != "" and name != "":
    dic[name] = tmp_contig

for k in sorted(dic, key=lambda k: len(dic[k]), reverse=True):
    print(k)
    print(textwrap.fill(dic[k], 100, break_on_hyphens=False))