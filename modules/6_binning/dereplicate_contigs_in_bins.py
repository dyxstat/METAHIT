#!/usr/bin/env python
import os
import sys


def load_scores(stats_path):
    scores = {}
    with open(stats_path) as handle:
        for line in handle:
            if "completeness" in line:
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 6:
                continue
            scores[cut[0]] = float(cut[1]) - 5 * float(cut[2]) + 0.0000000001 * float(cut[5])
    return scores


def bin_name_from_file(path):
    return ".".join(os.path.basename(path).split(".")[:-1])


stats_file, bins_folder, out_folder = sys.argv[1:4]
remove_all_duplicates = len(sys.argv) > 4 and sys.argv[4] == "remove"

print("Loading in bin completion and contamination scores...")
bin_scores = load_scores(stats_file)

print("Loading in contigs in each bin...")
contig_mapping = {}
for bin_file in os.listdir(bins_folder):
    bin_name = bin_name_from_file(bin_file)
    with open(os.path.join(bins_folder, bin_file)) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            contig = line[1:].strip()
            if contig not in contig_mapping:
                contig_mapping[contig] = bin_name
            elif remove_all_duplicates:
                contig_mapping[contig] = None
            elif bin_scores[bin_name] > bin_scores[contig_mapping[contig]]:
                contig_mapping[contig] = bin_name

print("Making a new dereplicated version of each bin file")
os.makedirs(out_folder, exist_ok=True)
for bin_file in os.listdir(bins_folder):
    bin_name = bin_name_from_file(bin_file)
    in_path = os.path.join(bins_folder, bin_file)
    out_path = os.path.join(out_folder, bin_file)
    at_least_one = False
    store = False

    with open(in_path) as src, open(out_path, "w") as out:
        for line in src:
            if line.startswith(">"):
                contig = line[1:].strip()
                store = contig_mapping.get(contig) == bin_name
                at_least_one = at_least_one or store
            if store:
                out.write(line)

    if not at_least_one:
        os.remove(out_path)
