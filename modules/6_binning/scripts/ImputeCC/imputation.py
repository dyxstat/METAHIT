#!/usr/bin/env python
# coding: utf-8

import logging
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as scisp
from .utility import random_walk_cpu, get_contigs_with_marker_genes



##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")
# package logger
logger = logging.getLogger(__name__)


'''
Cut off selection
'''
class ImputeMatrix:
    def __init__(self, contig_file, matrix_file, marker_file, gene_cov, rwr_rp, rwr_thres, max_markers):
        self.rwr_rp = rwr_rp #probability of restart
        self.rwr_perc = rwr_thres #threshold to remove the imputed Hi-C contacts to keep the sparsity of Hi-C matrix
   
        '''
        Information of assembled contigs
        dict_contig_len:    dict[contig_name] = contig_length
        dict_contigRev:     dict[contig_name] = row number in contig info np array
        '''
        ######Input the contig information
        #####input the blastn result for the contigs on reference genome###########

        self.contig_info = pd.read_csv(contig_file , sep = ',' , header = 0).values

        self.dict_contig_len = {}
        for i in self.contig_info:
            self.dict_contig_len[i[0]] = i[2]
            
        self.dict_contigRev = {} #dict[contig_name] = row number in contig info np array
        for i , j in enumerate(self.contig_info[: , 0]):
            self.dict_contigRev[j] = i
            
        ######Input normalized Hi-C matrix#######
        self.normcc_matrix = scisp.load_npz(matrix_file)
        self.normcc_matrix.tolil()
        
        '''
        Detect Marker genes
        '''
        ####Input marker gene information
        
            
        self.marker_contigs, self.marker_contig_counts, self.contig_markers = get_contigs_with_marker_genes(marker_file , gene_cov, self.dict_contigRev)
        if len(self.contig_markers) > max_markers:
            contig_markers_len = []
            contig_markers_name = []
            for contig in self.contig_markers.keys():
                contig_markers_name.append(contig)
                contig_markers_len.append(self.dict_contig_len[contig])
    
            contig_markers_len = np.array(contig_markers_len)
            contig_markers_name = np.array(contig_markers_name)
            
            retained_contigs = contig_markers_name[np.flip(np.argsort(contig_markers_len)[-max_markers:])]
            contig_markers_retained = {}
            
            for contig in retained_contigs:
                contig_markers_retained[contig] = self.contig_markers[contig]

            del self.marker_contigs, self.marker_contig_counts, self.contig_markers
            
            self.contig_markers = contig_markers_retained
            del contig_markers_retained

            self.marker_contigs = {}
            for i, j in self.contig_markers.items():
                for k in j:
                    if k not in self.marker_contigs:
                        self.marker_contigs[k] = [i]
                    else:
                        self.marker_contigs[k].append(i)
                        
            self.marker_contig_counts = {}
            for i, j in self.marker_contigs.items():
                self.marker_contig_counts[i] = len(j)
        logger.info('Conduct CRWR starting from {} marker-gene-containing contigs.'.format(len(self.contig_markers)))                      
        self.contig_local, self.dict_contigRevLocal, self.imputed_matrix = self._imputation()
        


    def  _imputation(self):
        '''
        Impute normcc-normalized Hi-C contact matrix using randow walk with restart
        '''
        ###Only impute part of matrix instead of the whole matrix
        _contig_with_marker = list(self.contig_markers.keys())
        _num_marker_contig = len(_contig_with_marker)
        
        _contig_with_marker_index = list(np.sort([self.dict_contigRev[i] for i in _contig_with_marker]))
        _contig_local = self.contig_info[_contig_with_marker_index, :]
        
        _dict_contigRevLocal = {}
        for i, j in enumerate(_contig_local[: , 0]):
            _dict_contigRevLocal[j] = i
        
        for i in range(self.contig_info.shape[0]):
            if i not in _contig_with_marker_index:
                _contig_with_marker_index.append(i)
                
        _normcc_shuffle = self.normcc_matrix.copy().tocsr()[_contig_with_marker_index, :]
        _normcc_shuffle = _normcc_shuffle.tocsc()[: , _contig_with_marker_index]

        A = _normcc_shuffle.copy()
        A = A - scisp.diags(A.diagonal())
        B = A + scisp.diags((A.sum(axis=0).A.ravel() == 0).astype(int))
        d = scisp.diags(1 / B.sum(axis=0).A.ravel())
        P = d.dot(B).astype(float)
        Q = random_walk_cpu(P , self.rwr_rp , self.rwr_perc , 0.01 , _num_marker_contig)

        E = Q.copy()
        E += E.T
        d = E.sum(axis=0).A.ravel()
        d[d == 0] = 1
        b = scisp.diags(1 / np.sqrt(d))
        E = b.dot(E).dot(b)
        
        return _contig_local, _dict_contigRevLocal, E.tolil()
        
        
    
    
