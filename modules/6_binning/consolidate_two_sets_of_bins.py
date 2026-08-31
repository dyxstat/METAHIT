#!/usr/bin/env python
import os
import shutil
import sys


def load_good_bins(stats_path, min_completion, max_contamination):
    good = {}
    with open(stats_path) as handle:
        for line in handle:
            if "completeness" in line:
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 3:
                continue
            if float(cut[1]) > min_completion and float(cut[2]) < max_contamination:
                good[cut[0] + ".fa"] = None
    return good


def load_stats(stats_path):
    stats = {}
    summary = {}
    with open(stats_path) as handle:
        for line in handle:
            if "completeness" in line:
                summary["header"] = line
                continue
            cut = line.rstrip("\n").split("\t")
            stats[cut[0] + ".fa"] = (float(cut[1]), float(cut[2]))
            summary[cut[0] + ".fa"] = line
    return stats, summary


def contig_lengths(bin_folder, bin_files):
    bins = {}
    for bin_file in bin_files:
        lengths = {}
        contig_len = 0
        contig_name = ""
        with open(os.path.join(bin_folder, bin_file)) as handle:
            for line in handle:
                if line.startswith(">"):
                    if contig_name:
                        lengths[contig_name] = contig_len
                    contig_name = line[1:].strip()
                    contig_len = 0
                else:
                    contig_len += len(line.strip())
        if contig_name:
            lengths[contig_name] = contig_len
        bins[bin_file] = lengths
    return bins


def overlap_ratio(lengths_a, lengths_b):
    match_a = 0
    match_b = 0
    mismatch_a = 0
    mismatch_b = 0

    for contig, length in lengths_a.items():
        if contig in lengths_b:
            match_a += lengths_b[contig]
        else:
            mismatch_a += length

    for contig, length in lengths_b.items():
        if contig in lengths_a:
            match_b += lengths_a[contig]
        else:
            mismatch_b += length

    denom_a = match_a + mismatch_a
    denom_b = match_b + mismatch_b
    ratio_a = 100 * match_a / denom_a if denom_a else 0
    ratio_b = 100 * match_b / denom_b if denom_b else 0
    return max(ratio_a, ratio_b)


bin_folder_1, bin_folder_2, stats_file_1, stats_file_2, output_folder = sys.argv[1:6]
min_completion = float(sys.argv[6])
max_contamination = float(sys.argv[7])
contamination_penalty = float(sys.argv[8]) if len(sys.argv) > 8 else 5.0

print(f"Loading list of good bins (comp>{min_completion}%, cont<{max_contamination}%)")
good_bins_1 = load_good_bins(stats_file_1, min_completion, max_contamination)
good_bins_2 = load_good_bins(stats_file_2, min_completion, max_contamination)

print("load in the info about the contigs in each bin...")
bins_1 = contig_lengths(bin_folder_1, good_bins_1)
bins_2 = contig_lengths(bin_folder_2, good_bins_2)

print("make all possible comparisons between the two bin sets, and record total % identical length")
all_bin_pairs = {}
for bin_1, lengths_1 in bins_1.items():
    all_bin_pairs[bin_1] = {}
    for bin_2, lengths_2 in bins_2.items():
        all_bin_pairs[bin_1][bin_2] = overlap_ratio(lengths_1, lengths_2)

print("load in completion and contamination scores of all the bins")
bins_1_stats, bins_1_summary = load_stats(stats_file_1)
bins_2_stats, bins_2_summary = load_stats(stats_file_2)

print("go through first group, pull out identical bins from second group, and choose best")
os.makedirs(output_folder, exist_ok=True)
new_summary_file = bins_1_summary["header"]
bins_2_matches = {}
bin_ct = 1

for bin_1, pairs in all_bin_pairs.items():
    score = bins_1_stats[bin_1][0] - bins_1_stats[bin_1][1] * contamination_penalty
    found_better = False
    for bin_2, overlap in pairs.items():
        if overlap < 80:
            continue
        bins_2_matches[bin_2] = None
        if bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1] * contamination_penalty > score:
            shutil.copy2(os.path.join(bin_folder_2, bin_2), os.path.join(output_folder, f"bin.{bin_ct}.fa"))
            new_summary_file += "bin." + str(bin_ct) + "\t" + "\t".join(bins_2_summary[bin_2].split("\t")[1:])
            found_better = True
    if not found_better:
        shutil.copy2(os.path.join(bin_folder_1, bin_1), os.path.join(output_folder, f"bin.{bin_ct}.fa"))
        new_summary_file += "bin." + str(bin_ct) + "\t" + "\t".join(bins_1_summary[bin_1].split("\t")[1:])
    bin_ct += 1

print("retrieve bins from second group that were not found in first group")
for bin_2, stats in bins_2_stats.items():
    if stats[0] < min_completion or stats[1] > max_contamination:
        continue
    if bin_2 in bins_2_matches:
        continue
    shutil.copy2(os.path.join(bin_folder_2, bin_2), os.path.join(output_folder, f"bin.{bin_ct}.fa"))
    new_summary_file += "bin." + str(bin_ct) + "\t" + "\t".join(bins_2_summary[bin_2].split("\t")[1:])
    bin_ct += 1

with open(output_folder + ".stats", "w") as out:
    out.write(new_summary_file)

print("There were " + str(bin_ct) + " bins cherry-picked from the original sets!")
