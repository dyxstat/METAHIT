#!/usr/bin/env python3

import os
import argparse
import gzip
import pickle
import logging
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import font_manager
import matplotlib.pyplot as plt
import seaborn
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.sparse import load_npz

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("main")
LABEL_FONT_SIZE = 42
FONT_FAMILY = "Arial"


def require_arial():
    root = Path(__file__).resolve().parents[2]
    font_dir = root / "conda_envs" / "metahict_env" / "fonts"
    arial_files = [font_dir / name for name in ("arial.ttf", "arialbd.ttf", "arialbi.ttf", "ariali.ttf")]
    missing = [path for path in arial_files if not path.is_file()]
    if missing:
        raise RuntimeError(f"Arial font file is missing: {missing[0]}")
    for path in arial_files:
        font_manager.fontManager.addfont(str(path))
    font_manager.findfont(FONT_FAMILY, fallback_to_default=False)
    matplotlib.rcParams.update(
        {
            "font.family": [FONT_FAMILY],
            "font.sans-serif": [FONT_FAMILY],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


require_arial()

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
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.25)
    heat = seaborn.heatmap(
        mat,
        square=True,
        cmap="rocket",
        ax=ax,
        cbar=True,
        cbar_ax=cax,
        linewidths=0,
    )
    cbar = heat.collections[0].colorbar
    cbar.ax.tick_params(labelsize=LABEL_FONT_SIZE)
    ax.hlines(tick_locs, *ax.get_xlim(), color='grey', linewidth=0.5, linestyle='-.')
    ax.vlines(tick_locs, *ax.get_ylim(), color='grey', linewidth=0.5, linestyle='-.')

    # Remove bin labels entirely, keep only boundary lines
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    plt.tight_layout(pad=1.5)
    out_png = os.path.join(outdir, "heatmap.png")
    out_pdf = os.path.join(outdir, "heatmap.pdf")
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    plt.savefig(out_pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    plt.close()
    logger.info(f"Saved heatmap to: {out_png}")
    logger.info(f"Saved heatmap to: {out_pdf}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contact-map", required=True, help="Path to .npz contact matrix")
    parser.add_argument("--contig-info", required=True, help="MetaCC contig information CSV")
    parser.add_argument("--BIN", required=True, help="Path to final_bins.p.gz")
    parser.add_argument("--OUTDIR", required=True, help="Output directory")
    parser.add_argument("--max_image", type=int, default=5000, help="Max image size (default: 5000)")
    args = parser.parse_args()

    logger.info("Loading sparse matrix...")
    sparse_mat = load_npz(args.contact_map)

    logger.info("Loading bin mapping...")
    with gzip.open(args.BIN, 'rb') as f:
        contig_to_bin = pickle.load(f)

    logger.info("Loading MetaCC contig order...")
    contig_info = pd.read_csv(args.contig_info)
    if "name" not in contig_info.columns:
        raise ValueError(f"Missing 'name' column in {args.contig_info}")
    if sparse_mat.shape[0] != len(contig_info):
        raise ValueError(
            f"Contact matrix dimension {sparse_mat.shape[0]} does not match "
            f"the {len(contig_info)} MetaCC contig rows"
        )
    order = range(len(contig_info))
    idx_to_name = dict(enumerate(contig_info["name"].astype(str)))

    plot(sparse_mat, contig_to_bin, order, idx_to_name, args.OUTDIR, max_image_size=args.max_image)

if __name__ == "__main__":
    main()
