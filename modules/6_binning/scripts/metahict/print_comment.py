#!/usr/bin/env python
import sys

comm = sys.argv[1]
delim = sys.argv[2]
width = 120
max_len = 90

print("\n" + delim * width)

line = ""
for word in comm.split():
    if len(line) + 1 + len(word) > max_len:
        edge1 = int((width - len(line)) / 2 - 5)
        edge2 = int(width - edge1 - len(line) - 10)
        print(delim * 5 + " " * edge1 + line + " " * edge2 + delim * 5)
        line = word
    else:
        line = (line + " " + word).strip()

edge1 = int((width - len(line)) / 2 - 5)
edge2 = int(width - edge1 - len(line) - 10)
print(delim * 5 + " " * edge1 + line + " " * edge2 + delim * 5)
print(delim * width + "\n")
