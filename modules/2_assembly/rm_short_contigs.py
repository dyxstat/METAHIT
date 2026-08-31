#!/usr/bin/env python
from __future__ import print_function

import argparse
import re


def fasta_records(path):
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def header_length(header):
    match = re.search(r"(?:^|_)length[_=](\d+)(?:_|$)", header)
    if match:
        return int(match.group(1))
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA records shorter than a minimum length."
    )
    parser.add_argument("min_len", type=int, help="Minimum contig length to retain")
    parser.add_argument("fasta", help="Input FASTA")
    args = parser.parse_args()

    for header, seq in fasta_records(args.fasta):
        length = header_length(header)
        if length is None:
            length = len(seq)
        if length < args.min_len:
            continue
        print(header)
        for i in range(0, len(seq), 80):
            print(seq[i:i + 80])


if __name__ == "__main__":
    main()
