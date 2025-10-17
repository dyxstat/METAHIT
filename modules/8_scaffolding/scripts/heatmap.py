#!/usr/bin/env python3

import os
import argparse
import gzip
import pickle
import logging
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import seaborn
from matplotlib import ticker
from scipy.sparse import load_npz

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("main")

def downsample(mat, factor):
    if factor <= 1:
        return mat
    size = mat.shape[0]
    new_size = int(np.ceil(size / factor))
    reduced = np.zeros((new_size, new_size), dtype=mat.dtype)
    for i in range(new_size):
        for j in range(new_size):
            block = mat[i*factor:(i+1)*factor, j*factor:(j+1)*factor]
            reduced[i, j] = block.mean() if block.size else 0
    return reduced

def plot(sparse_mat, contig_to_bin, order, idx_to_name, outdir, max_image_size=5000, width=25, height=22, dpi=300):
    os.makedirs(outdir, exist_ok=True)
    ordered_names = [idx_to_name[i] for i in order if i in idx_to_name]

    bins = defaultdict(list)
    for i, name in enumerate(ordered_names):
        if name in contig_to_bin:
            bins[contig_to_bin[name]].append(i)

    bin_ids = sorted(bins)
    tick_locs, final_indices = [0], []
    for b in bin_ids:
        idxs = bins[b]
        final_indices.extend(idxs)
        tick_locs.append(tick_locs[-1] + len(idxs))

    logger.info(f"Total contigs in plot: {len(final_indices)}")

    # Extract only the needed submatrix and densify that
    mat = sparse_mat[final_indices, :][:, final_indices].toarray()

    if max(mat.shape) > max_image_size:
        factor = int(np.ceil(max(mat.shape) / float(max_image_size)))
        logger.info(f"Downsampling matrix from {mat.shape} by factor {factor}")
        mat = downsample(mat, factor)
        tick_locs = np.floor(np.array(tick_locs) / factor).astype(int)

    np.fill_diagonal(mat, 0)
    mat = np.log(mat + 0.01)

    fig, ax = plt.subplots(figsize=(width, height))
    seaborn.heatmap(mat, square=True, cmap="rocket", ax=ax, cbar=True, linewidths=0)
    ax.hlines(tick_locs, *ax.get_xlim(), color='grey', linewidth=0.5, linestyle='-.')
    ax.vlines(tick_locs, *ax.get_ylim(), color='grey', linewidth=0.5, linestyle='-.')

    # Remove bin labels entirely, keep only boundary lines
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    plt.tight_layout()
    out_path = os.path.join(outdir, "heatmap.png")
    plt.savefig(out_path, dpi=dpi)
    plt.close()
    logger.info(f"Saved heatmap to: {out_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contact-map", required=True, help="Path to .npz contact matrix")
    parser.add_argument("--BIN", required=True, help="Path to final_bins.p.gz")
    parser.add_argument("--ORDER", required=True, help="Path to contact_map.p.gz")
    parser.add_argument("--OUTDIR", required=True, help="Output directory")
    parser.add_argument("--max_image", type=int, default=5000, help="Max image size (default: 5000)")
    args = parser.parse_args()

    logger.info("Loading sparse matrix...")
    sparse_mat = load_npz(args.contact_map)

    logger.info("Loading bin mapping...")
    with gzip.open(args.BIN, 'rb') as f:
        contig_to_bin = pickle.load(f)

    logger.info("Loading order and seq info...")
    with gzip.open(args.ORDER, 'rb') as f:
        obj = pickle.load(f)
        order = obj.order.order["pos"]
        idx_to_name = {s.localid_metacc: s.name for s in obj.seq_info_metacc}

    plot(sparse_mat, contig_to_bin, order, idx_to_name, args.OUTDIR, max_image_size=args.max_image)

if __name__ == "__main__":
    main()
