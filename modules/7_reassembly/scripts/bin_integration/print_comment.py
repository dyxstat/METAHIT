#!/usr/bin/env python
# This script prints any comment in a structured and prety way.
from __future__ import print_function
import sys
comm = sys.argv[1]
delim = sys.argv[2]

print('\n'+delim*120)

max_len = 90

cut = comm.split(" ")
line = ""
for word in cut:
    if (len(line) + 1 + len(word)) > max_len:
        edge1 = (120-len(line))//2 - 5
        edge2 = 120-edge1-len(line) - 10
        print(delim*int(5) + " "*int(edge1) +
              line + " "*int(edge2) + delim*int(5))
        line = word
    else:
        line = line+" "+word
edge1 = (120-len(line))//2 - 5
edge2 = 120-edge1-len(line) - 10
print(delim*int(5) + " "*int(edge1) + line + " "*int(edge2) + delim*int(5))

print(delim*120+'\n')
