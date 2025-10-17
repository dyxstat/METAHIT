#!/usr/bin/env python
# coding: utf-8
import warnings
##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")
import logging
import numpy as np
import pandas as pd
import leidenalg
import operator
import scipy.sparse as scisp
import igraph as ig
import leidenalg
import os
from multiprocessing import Pool

# External script
from .pre_markers import Markers
from .utility import match_contigs
from .ensemble import Ensemble


##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")
# package logger
logger = logging.getLogger(__name__)


'''
Cut off selection
'''
class FinalCluster:
    def __init__(self, contig_info, contig_local, 
                 dict_contigRev, dict_contigRevLocal, dict_contig_len, 
                 contig_markers, bins, bin_of_contigs,
                 normcc_matrix, imputed_matrix, 
                 bac_mg_table, ar_mg_table, path, 
                 intra, inter, cont_weight, 
                 min_comp, max_cont, report_quality, min_binsize, n_process=12):
        
        self.contig_info = contig_info
        self.contig_local = contig_local
        self.dict_contigRev = dict_contigRev
        self.dict_contigRevLocal = dict_contigRevLocal
        self.dict_contig_len = dict_contig_len
        ### Results from get_contigs_with_marker_genes in marker file
        self.contig_markers = contig_markers
        
        self.normcc_matrix = normcc_matrix
        self.imputed_matrix = imputed_matrix

        self.w_intra = intra
        self.w_inter = inter
        self.path = path
        self.minbinsize = min_binsize
        
        self.markers = Markers()
        _ = self.markers.marker_gene_tables(bac_mg_table, ar_mg_table) 
        
        self.bins = bins 
        self.bin_of_contigs = bin_of_contigs
        
        #Construct the input graph to Leiden clistering
        map_matrix = self.normcc_matrix.copy().tocoo()
        vcount = map_matrix.shape[0]
        sources = map_matrix.row
        targets = map_matrix.col
        wei = map_matrix.data
        index = sources<targets
        sources = sources[index]
        targets = targets[index]
        wei = wei[index]
        edgelist = list(zip(sources, targets))
        g = ig.Graph(vcount, edgelist)
        del map_matrix, vcount, sources, targets, index, edgelist 
        
        res_option = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        
        logger.info('Begin Leiden clustering with %d processes...' % (n_process))

        #CLustering without using information of marker genes
        p = Pool(n_process)
        total_5_nofix = p.starmap(self.leiden_clustering_without_marker_genes, [(g, leidenalg.RBConfigurationVertexPartition, wei, res) for res in res_option])
        p.close()
        
        ##Clustering without using marker genes
        res_para_nofix = res_option[total_5_nofix.index(max(total_5_nofix))]
        part_nofix = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=wei , resolution_parameter = res_para_nofix , n_iterations = -1)
        part_nofix = list(part_nofix)

        bin2contig_nofix_ori = {}
        for ci in range(len(part_nofix)):
            if sum(self.contig_info[part_nofix[ci] , 2]) >= self.minbinsize:
                for id in part_nofix[ci]:
                    if 'BIN'+ str(ci) in bin2contig_nofix_ori:
                        bin2contig_nofix_ori['BIN'+str(ci)].append(self.contig_info[id , 0])
                    else:
                        bin2contig_nofix_ori['BIN'+str(ci)] = [self.contig_info[id , 0]]
        
        bin2contig_nofix_decontamin_first_time = self.decontamin(bin2contig_nofix_ori)
        self.bin2contig_nofix_decontamin = self.decontamin(bin2contig_nofix_decontamin_first_time)
        
        #num_50 = 0
        #for i in self.bin2contig_nofix_decontamin.values():
        #    _, comp, cont = self.markers.bin_quality(i)
        #    if comp >= 50 and cont <= 5:
        #        num_50 += 1
        #logger.info('There are {} high quality MAGs for nofix'.format(num_50))
        
        ##Clustering using marker genes
        '''
        Semi-supervised Leiden
        '''
        is_membership_fixed = []
        new_membership = []
        index = len(self.bins)

        for i in self.contig_info[: , 0]:
            if i in self.bin_of_contigs:
                new_membership.append(self.bin_of_contigs[i])
                is_membership_fixed.append(True)
            else:
                new_membership.append(index)
                index += 1
                is_membership_fixed.append(False)

        p = Pool(n_process)
        total_5_fix = p.starmap(self.leiden_clustering_with_marker_genes, [(g, new_membership, wei, res, is_membership_fixed) for res in res_option])
        p.close()
        
        ##Clustering without using marker genes
        res_para_fix = res_option[total_5_fix.index(max(total_5_fix))]

        new_partition = leidenalg.RBConfigurationVertexPartition(g, new_membership, weights=wei , resolution_parameter = res_para_fix)
        optimiser = leidenalg.Optimiser()
        _ = optimiser.optimise_partition(new_partition, is_membership_fixed=is_membership_fixed, n_iterations = -1)
        mem = new_partition.membership

        part_fix = [[] for i in range(len(set(mem)))]
        for i, j in enumerate(mem):
            part_fix[j].append(i)

        bin2contig_fix_ori = {}
        for ci in range(len(part_fix)):
            if sum(self.contig_info[part_fix[ci] , 2]) >= self.minbinsize:
                for id in part_fix[ci]:
                    if 'BIN' + str(ci) in bin2contig_fix_ori :
                        bin2contig_fix_ori ['BIN' + str(ci)].append(self.contig_info[id , 0])
                    else:
                        bin2contig_fix_ori ['BIN' + str(ci)] = [self.contig_info[id , 0]]
        
        #Decontaminate for only one time
        #self.bin2contig_fix_decontamin = self.decontamin(bin2contig_fix_ori)
        
        #Two times decontamination
        bin2contig_fix_decontamin_first_time = self.decontamin(bin2contig_fix_ori)
        self.bin2contig_fix_decontamin = self.decontamin(bin2contig_fix_decontamin_first_time)
        logger.info('Leiden clustering finished.')
        
        #num_50 = 0
        #for i in self.bin2contig_fix_decontamin.values():
        #    _, comp, cont = self.markers.bin_quality(i)
        #    if comp >= 50 and cont <= 5:
        #        num_50 += 1
        #logger.info('There are {} high quality MAGs for prefix'.format(num_50))
        
        ## Ensemble step #######
        """
        d[binning method][bin ID] -> set(cid1, cid2, ... cidN)
        Contigs for bins across all binning methods.
        """
        # transfer the current format to what the ensemble method needs
        bins_combine = {}
        bins_combine['fix'] = self.bin2contig_fix_decontamin.copy()
        bins_combine['nofix'] = self.bin2contig_nofix_decontamin.copy()
        
        methods_sorted = sorted(bins_combine.keys())
        
        ## Ensemble parameter
        em_mode = 'greedy' 
        #no_bin_matching = True  # greedy
        
        CONT_WEIGHT = cont_weight  # help="weight given to contamination for assessing genome quality", type=float, default=2)
        SEL_MIN_COMP = min_comp  # add_argument('-x', '--sel_min_comp', help="minimum completeness of bin to consider during bin selection process", type=float, default=50)
        SEL_MAX_CONT = max_cont  # add_argument('-y', '--sel_max_cont', help="maximum contamination of bin to consider during bin selection process", type=float, default=10)
        REPORT_MIN_QUALITY = report_quality  # add_argument('--report_min_quality', help="minimum quality of bin to report", type=float, default=10)

        BIN_PREFIX = 'bin'
        ensemble = Ensemble()
        # ensemble of binning results

        bin_num = 0
        total_comp = 0
        total_cont = 0
        total_quality = 0
        out_bin_quality = {}
        
        ######################################################################
        ########################## greedy  ##############################
        ######################################################################
        self.bin2contig_ensemble = {}

        while True:
            bins_quality = ensemble._bin_quality(self.markers, bins_combine, self.dict_contig_len, methods_sorted, CONT_WEIGHT)
            matched_sets = ensemble._matched_bin_sets(bins_combine, bins_quality)
            if len(matched_sets) == 0:
                break  # no bins to be resolve
            
            new_bin, _, _ = ensemble._reconstruct_match_sets(matched_sets[0],
                                                            bins_combine,
                                                            self.dict_contig_len,
                                                            em_mode=em_mode)
            

            _, comp, cont = self.markers.bin_quality(new_bin.keys())
            cc_score = comp - cont * CONT_WEIGHT

            if cc_score < REPORT_MIN_QUALITY:
                break
            else:
                if (comp >= SEL_MIN_COMP) and (cont <= SEL_MAX_CONT):
                    total_comp += comp
                    total_cont += cont
                    total_quality += cc_score

                    # report selection
                    bin_num += 1
                    out_bin_in = '%s_%d' % (BIN_PREFIX, bin_num)
                    out_bin_quality[out_bin_in] = (comp, cont)
                    primary_bm, primary_bid, _q, _n50, _gs = matched_sets[0][0]

                    # write out contig info
                    matched_bins = {}
                    for m in matched_sets[0]:
                        matched_bins[m[0]] = m[1]
                    
                    self.bin2contig_ensemble[bin_num] = list(new_bin.keys())
                        
                ensemble._update_bins(bins_combine, new_bin.keys())
                
        num_50 = 0
        for i in self.bin2contig_ensemble.values():
            _, comp, cont = self.markers.bin_quality(i)
            if comp >= 50 and cont <= 5:
                num_50 += 1
        #logger.info('After the ensemble, there are {} high quality MAGs'.format(num_50))
                
                    
        with open(os.path.join(self.path , 'tmp' , 'cluster_imputecc.txt'),'w') as out:
            for key , value in self.bin2contig_ensemble.items():
                for contig in value:
                    out.write(str(contig)+ '\t' +str(key)+ '\n')           
      
    def leiden_clustering_without_marker_genes(self, g, partition, wei, res):
        # part_nofix = list(leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=wei , resolution_parameter = res , n_iterations = -1))
        part_nofix = list(leidenalg.find_partition(g , partition , weights=wei , resolution_parameter = res , n_iterations = -1))
        cluster_nofix = []
        for ci in range(len(part_nofix)):
            if sum(self.contig_info[part_nofix[ci] , 2]) >= self.minbinsize:
                temp = []
                for id in part_nofix[ci]:
                    temp.append(self.contig_info[id , 0])
                cluster_nofix.append(temp)
 
        num_50 = 0
        for i in cluster_nofix:
            _, comp, cont = self.markers.bin_quality(i)
            if cont <= 5:
                if comp >= 90:
                    num_50 += 3
                elif comp >= 70:
                    num_50 += 2
                elif comp >= 50:
                    num_50 += 1
        return num_50

    def leiden_clustering_with_marker_genes(self, g, membership, wei, res, is_membership_fixed):
        new_partition = leidenalg.RBConfigurationVertexPartition(g, membership , weights=wei , resolution_parameter = res)
        optimiser = leidenalg.Optimiser()
        _ = optimiser.optimise_partition(new_partition, is_membership_fixed=is_membership_fixed, n_iterations = -1)
        mem = new_partition.membership
        
        part_fix = [[] for i in range(len(set(mem)))]
        for i, j in enumerate(mem):
            part_fix[j].append(i)

        cluster_fix = []
        for ci in range(len(part_fix)):
            if sum(self.contig_info[part_fix[ci] , 2]) >= self.minbinsize:
                temp = []
                for id in part_fix[ci]:
                    temp.append(self.contig_info[id , 0])
                cluster_fix.append(temp)
 
        num_50 = 0
        for i in cluster_fix:
            _, comp, cont = self.markers.bin_quality(i)
            if cont <= 5:
                if comp >= 90:
                    num_50 += 3
                elif comp >= 70:
                    num_50 += 2
                elif comp >= 50:
                    num_50 += 1
        return num_50
    
    def decontamin(self, bin2contig_ori):
        _contam_bins = []
        _decontam_bin = {}
        _bin_num_decontam = 0
        for key, value in bin2contig_ori.items():
            _, _comp, _cont = self.markers.bin_quality(value)
            if _comp >= 50 and _cont > 10:
                _contam_bins.append(key)
            else:
                _bin_name = 'BIN' + str(_bin_num_decontam)
                _decontam_bin[_bin_name] = value
                _bin_num_decontam += 1
        
        if len(_contam_bins) > 0:                
            for _contam_bin in _contam_bins:
                _contig_in_contam = bin2contig_ori[_contam_bin]  # find all contigs in contamined bins

                _index_marker_contig_local = []
                _contig_markers_local = {}
                
                for _contig in _contig_in_contam:
                    if _contig in self.dict_contigRevLocal: #dict_contigRevLocal shows the reverse index of contigs in imputed matrix
                        _index_marker_contig_local.append(self.dict_contigRevLocal[_contig])
                        _contig_markers_local[_contig] = self.contig_markers[_contig] 
                    
                if len(_index_marker_contig_local) > 0:
                    #Find imputed matrix corresponding to the contaminated bins
                    _map_marker_contig_contam = self.imputed_matrix.copy().tocsr()[_index_marker_contig_local, :]
                    _map_marker_contig_contam = _map_marker_contig_contam.tocsc()[: , _index_marker_contig_local]

                    _contig_info_marker_contig_contam = self.contig_local[_index_marker_contig_local, : ]
                    _dict_contigRevContam = {}
                    for i, j in enumerate(_contig_info_marker_contig_contam[: , 0]):
                        _dict_contigRevContam[j] = i
                        
                    #Change the format of _contig_markers_local
                    _marker_contigs_local = {}
                    _marker_contig_counts_local = {}
                                        
                    for i , j in _contig_markers_local.items():
                        for k in j:
                            if k in _marker_contigs_local:
                                _marker_contigs_local[k].append(i)
                            else:
                                _marker_contigs_local[k] = [i]
                                
                    for i, j in _marker_contigs_local.items():
                        _marker_contig_counts_local[i] = len(j)
                        
                    _my_gene_counts_local = list(_marker_contig_counts_local.values())
                    _my_gene_counts_local.sort(reverse=True)

                    _smg_iteration = {}

                    n = 0
                    _unique_my_gene_counts = list(set(_my_gene_counts_local))
                    _unique_my_gene_counts.sort(reverse=True)

                    # Get contigs for each iteration of single-copy marker gene
                    for _g_count in _unique_my_gene_counts:
                        # Get the single-copy marker genes with maximum count of contigs and
                        # sort them in the descending order of the total marker genes contained

                        _total_contig_mgs = {}
                        
                        ####item here is the single copy marker genes
                        #If multiple genes correspond to the same number of contigs, we will sort the genes according to the number of genes on contigs
                        for _item in _marker_contig_counts_local:
                            if _marker_contig_counts_local[_item] == _g_count:
                                _total_contig_lengths = 0
                                
                                #Marker contigs are dict gene: contigs###
                                for _contig in _marker_contigs_local[_item]:
                                    _contig_mg_counts = len(_contig_markers_local[_contig])
                                    _total_contig_lengths += _contig_mg_counts

                                _total_contig_mgs[_item] = _total_contig_lengths

                        _total_contig_mgs_sorted = sorted(
                            _total_contig_mgs.items(), key=operator.itemgetter(1), reverse=True
                        )

                        for _item in _total_contig_mgs_sorted:
                            _smg_iteration[n] = _marker_contigs_local[_item[0]]
                            n += 1
                            
                    if len(self.imputed_matrix.tocoo().data) == 0:
                        _bins_contam = {}
                        _bin_of_contigs_contam = {}
                        for i in range(len(_smg_iteration[0])):
                            contig_num = _smg_iteration[0][i]
                            _bins_contam[i] = [contig_num]
                            _bin_of_contigs_contam[contig_num] = i
                    else:    
                        _bins_contam, _bin_of_contigs_contam, _, _, _ = match_contigs(
                        _smg_iteration,
                        _contig_markers_local,
                        _map_marker_contig_contam.tolil(),
                        _dict_contigRevContam,
                        np.percentile(self.imputed_matrix.tocoo().data , self.w_intra),
                        np.percentile(self.imputed_matrix.tocoo().data , self.w_inter)
                        )
                else:
                    _bin_of_contigs_contam = {}
                    _bins_contam = {}
                    

                #index_all_contigs_contam is for all contigs in the contaminated bins, not only the contigs with marker genes
                _index_all_contigs_contam = []
                for _contig in _contig_in_contam:
                    if _contig in self.dict_contigRev:
                        _index_all_contigs_contam.append(self.dict_contigRev[_contig])
                        
                _map_all_contigs_contam = self.normcc_matrix.copy().tocsr()[_index_all_contigs_contam , :]
                _map_all_contigs_contam = _map_all_contigs_contam.tocsc()[:, _index_all_contigs_contam]

                _contig_info_all_contam = self.contig_info[_index_all_contigs_contam, :]
                
                _is_membership_fixed = []
                _new_membership = []
                _index_membership = len(_bins_contam)
                for i in _contig_info_all_contam[: , 0]:
                    if i in _bin_of_contigs_contam:
                        _new_membership.append(_bin_of_contigs_contam[i])
                        _is_membership_fixed.append(True)
                    else:
                        _new_membership.append(_index_membership)
                        _index_membership += 1
                        _is_membership_fixed.append(False)   

                            
                                               
                _map_del = _map_all_contigs_contam.copy().tocoo()
                _vcount = _map_del.shape[0]
                _sources = _map_del.row
                _targets = _map_del.col
                _wei = _map_del.data
                _index = _sources < _targets
                _sources = _sources[_index]
                _targets = _targets[_index]
                _wei = _wei[_index]
                _edgelist = list(zip(_sources, _targets))
                _g = ig.Graph(_vcount, _edgelist)
                del _map_del, _vcount, _sources, _targets, _index, _edgelist

                _new_partition = leidenalg.RBConfigurationVertexPartition(_g, _new_membership, weights=_wei , resolution_parameter = 1)
                _optimiser = leidenalg.Optimiser()
                _optimiser.consider_empty_community = False
                _ = _optimiser.optimise_partition(_new_partition, is_membership_fixed=_is_membership_fixed, n_iterations = -1)
                _mem = _new_partition.membership

                _split_bin = [[] for i in range(len(set(_mem)))]
                for i, j in enumerate(_mem):
                    _split_bin[j].append(_contig_info_all_contam[i , 0])
        
                for i, j in enumerate(_split_bin):
                    gs = 0
                    for k in j:
                        gs += self.dict_contig_len[k]
                    if gs >= self.minbinsize:
                        _bin_name = 'BIN' + str(_bin_num_decontam)
                        _decontam_bin[_bin_name] = j
                        _bin_num_decontam += 1
        
        return _decontam_bin
