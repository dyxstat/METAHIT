import igraph as ig
import leidenalg
import os
import scipy.sparse as scisp
import pandas as pd
import numpy as np
import operator
from .utility import gen_bins, match_contigs
import logging
def PreCluster(marker_contig_counts, marker_contigs, contig_markers, imputed_matrix, dict_contigRevLocal, intra, inter):
    """
    Pre-clusters contigs based on marker gene counts and normalized contact matrix.

    Parameters:
        marker_contig_counts (dict): Counts of marker genes per contig.
        marker_contigs (dict): Mapping from marker genes to contigs.
        contig_markers (dict): Mapping from contigs to their marker genes.
        imputed_matrix (scipy.sparse matrix): Normalized contact matrix.
        dict_contigRevLocal (dict): Reverse mapping of contigs (local to global or vice versa).
        intra (float): Intra-bin percentile threshold.
        inter (float): Inter-bin percentile threshold.

    Returns:
        tuple: Bins dictionary and bin_of_contigs dictionary.
    """
    _my_gene_counts = list(marker_contig_counts.values())
    _my_gene_counts.sort(reverse=True)

    smg_iteration = {}
    n = 0
    _unique_my_gene_counts = sorted(list(set(_my_gene_counts)), reverse=True)

    # Get contigs for each iteration of single-copy marker gene
    for _g_count in _unique_my_gene_counts:
        total_contig_mgs = {}

        # Items here are the single-copy marker genes
        for item in marker_contig_counts:
            if marker_contig_counts[item] == _g_count:
                total_contig_lengths = 0

                # Marker contigs are dict gene: contigs
                for contig in marker_contigs[item]:
                    contig_mg_counts = len(contig_markers.get(contig, []))
                    total_contig_lengths += contig_mg_counts

                total_contig_mgs[item] = total_contig_lengths

        # Sort genes based on the number of marker genes they contain
        total_contig_mgs_sorted = sorted(
            total_contig_mgs.items(), key=operator.itemgetter(1), reverse=True
        )

        for item in total_contig_mgs_sorted:
            smg_iteration[n] = marker_contigs[item[0]]
            n += 1

    # Check if any marker genes were detected
    if not smg_iteration:
        logging.warning("No marker genes detected. Skipping pre-clustering.")
        return {}, {}

    # Proceed with clustering based on the presence of data in the imputed matrix
    if len(imputed_matrix.tocoo().data) == 0:
        bins = {}
        bin_of_contig = {}
        logging.info(f"Pre-clustering with {len(smg_iteration)} marker gene iterations.")

        for i in range(len(smg_iteration[0])):
            contig_num = smg_iteration[0][i]
            bins[i] = [contig_num]
            bin_of_contig[contig_num] = i
    else:
        bins, bin_of_contig, _, _, _ = match_contigs(
            smg_iteration,
            contig_markers,
            imputed_matrix.tolil(),
            dict_contigRevLocal,
            np.percentile(imputed_matrix.tocoo().data, intra),
            np.percentile(imputed_matrix.tocoo().data, inter)
        )

    logging.info(f"Pre-clustering finished with {len(bins)} preliminary bins established.")
    return bins, bin_of_contig


def Clust4CheckM(fasta_file, contig_info_file, normcc_matrix_file, path):
    _map_del = scisp.load_npz(normcc_matrix_file).tocoo()
    contig_info = pd.read_csv(contig_info_file , sep = ',' , header = 0).values
    
    #########Use Leiden Algorithm to do clustering########
    _vcount = _map_del.shape[0]
    _sources = _map_del.row
    _targets = _map_del.col
    _wei = _map_del.data
    _index = _sources<_targets
    _sources = _sources[_index]
    _targets = _targets[_index]
    _wei = _wei[_index]
    _edgelist = list(zip(_sources, _targets))
    g = ig.Graph(_vcount, _edgelist)
    
    #############determine the resolution parameter###########
    part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=_wei)
    part = list(part)
    _dict_cluster = {}
    # dict of communities

    for ci in range(len(part)):
        for id in part[ci]:
            _dict_cluster[contig_info[id , 0]] = 'group'+str(ci)
            
    with open(os.path.join(path , 'tmp' , 'cluster4checkm.txt'),'w') as out:
        for key , value in _dict_cluster.items():
            out.write(str(key)+ '\t' +str(value)+ '\n')
            
    with open(os.path.join(path , 'tmp' , 'dir4checkm.tsv'),'w') as out:
        out.write('INITIAL_BIN' + '\t' + os.path.join(path , 'tmp' ,'BIN4checkm'))
    
    gen_bins(fasta_file , os.path.join(path , 'tmp' , 'cluster4checkm.txt') , os.path.join(path , 'tmp' ,'BIN4checkm'))
    
    
            
