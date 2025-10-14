#!/usr/bin/env python
# coding: utf-8

'''
libraries for metacc
'''
from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from Bio.Restriction import Restriction
import Bio.SeqIO as SeqIO
import numpy as np
import pysam
import scipy.sparse as scisp
import tqdm
import os
from MetaCC.Script.utils import count_fasta_sequences, open_input
import logging

'''
libraries for bin3c
'''
import matplotlib
matplotlib.use('Agg')

from bin3C_python3.mzd import io_utils, sparse_utils
from bin3C_python3.mzd.seq_utils import *
from collections import OrderedDict, namedtuple
from numba import jit, int64, float64, void
import Bio.SeqIO as SeqIO
import logging
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pysam
import scipy.sparse as sp
import tqdm

# package logger
logger = logging.getLogger(__name__)

#######offset is the enumeration of the length######################
#localid is the index of the list#
#refid is the index of the fasta file, which is a global index# 

def mean_selector(name):
    """
    Basic Mean functions
    """
    def geometric_mean(x, y):
        return (x*y)**0.5

    def harmonic_mean(x, y):
        return 2*x*y/(x+y)

    def arithmetic_mean(x, y):
        return 0.5*(x+y)

    try:
        mean_switcher = {
            'geometric': geometric_mean,
            'harmonic': harmonic_mean,
            'arithmetic': arithmetic_mean
        }
        return mean_switcher[name]
    except KeyError:
        raise RuntimeError('unsupported mean type [{}]'.format(name))




SeqInfo_metacc = namedtuple('SeqInfo_metacc', ['localid_metacc', 'refid', 'name', 'sites', 'length', 'covcc']) #create a new class of tuple: SeqInfo

SeqInfo_bin3c = namedtuple('SeqInfo_bin3c', ['offset', 'refid', 'name', 'length', 'sites'])

class SiteCounter(object):

    def __init__(self, enzyme_names,  is_linear=True):
        """
        Simple class to count the total number of enzymatic cut sites for the given
        list if enzymes.
        :param enzyme_names: a list of enzyme names (proper case sensitive spelling a la NEB)
        :param is_linear: Treat sequence as linear.
        """
        if isinstance(enzyme_names, str):
            enzyme_names = [enzyme_names]
        assert (isinstance(enzyme_names, Iterable) and
                not isinstance(enzyme_names, str)), 'enzyme_names must of a collection of names'
        self.enzymes = [getattr(Restriction , en) for en in enzyme_names] ##get the enzyme name from the standard module
        self.is_linear = is_linear


    def count_sites(self, seq):
        """
        Count the number of sites found in the given sequence, where sites from
        all specified enzymes are combined

        :param seq: Bio.Seq object
        :return: the total number of sites
        """

        return sum(len(en.search(seq, self.is_linear)) for en in self.enzymes)




class Sparse2DAccumulator(object):
######create a 2D coo sparse matrix###########
    def __init__(self, N):
        self.shape = (N, N)
        self.mat = {}
        ##mat is a dictionary here
        # fixed counting type
        self.dtype = np.uint32

    def setitem(self, index, value):
        assert len(index) == 2 and index[0] >= 0 and index[1] >= 0, 'invalid index: {}'.format(index)
        self.mat[index] = value

    def getitem(self, index):
        if index in self.mat:
            return self.mat[index]
        else:
            return 0

    def get_coo(self, symm=True):
        """
        Create a COO format sparse representation of the accumulated values.

        :param symm: ensure matrix is symmetric on return
        :return: a scipy.coo_matrix sparse matrix
        """
        coords = [[], []]
        data = []
        m = self.mat
        for i, j in m.keys(): ##m.keys() will return a tuple of two values
            coords[0].append(i)
            coords[1].append(j)
            data.append(m[i, j])

        m = scisp.coo_matrix((data, coords), shape=self.shape, dtype=self.dtype)

        if symm:
            m += scisp.tril(m.T, k=-1)

        return m.tocoo()


def fast_norm_fullseq_bysite(rows, cols, data, sites):
    """
    In-place normalisation of the scipy.coo_matrix for full sequences

    :param rows: the COO matrix coordinate member variable (4xN array)
    :param cols: the COO matrix coordinate member variable (4xN array)
    :param data:  the COO matrix data member variable (1xN array)
    :param sites: per-element min(sequence_length, tip_size)
    """
    for n in range(data.shape[0]):
        i = rows[n]
        j = cols[n]
        data[n] *= 1.0/(sites[i] * sites[j])
        
class SeqOrder:

    FORWARD = 1
    REVERSE = -1

    ACCEPTED = True
    EXCLUDED = False

    STRUCT_TYPE = np.dtype([('pos', int), ('ori', int), ('mask', bool), ('length', int)])
    INDEX_TYPE = np.dtype([('index', int), ('ori', int)])

    def __init__(self, seq_info):
        """
        Initial order is determined by the order of supplied sequence information dictionary. Sequences
        are given surrogate ids using consecutive integers. Member functions expect surrogate ids
        not original names.

        The class also retains orientation and masking state. Orientation defines whether a sequence
        should be in its original direction (as read in) (1) or reverse complemented (-1).

        Masking state defines whether a input sequence shall been excluded from further consideration.
        (accepted=1, excluded=0)

        :param seq_info: sequence information dictionary
        """
        self._positions = None
        _ord = np.arange(len(seq_info), dtype=int)
        self.order = np.array(
            [(_ord[i], SeqOrder.FORWARD, SeqOrder.ACCEPTED, seq_info[i].length) for i in range(len(_ord))],
            dtype=SeqOrder.STRUCT_TYPE)

        self._update_positions()

    @staticmethod
    def asindex(_ord):
        """
        Convert a simple list or ndarray of indices, to INDEX_TYPE array with default forward orientation.

        :param _ord: list/ndarray of indices
        :return: INDEX_TYPE array
        """
        assert isinstance(_ord, (list, np.ndarray)), 'input must be a list or ndarray'
        return np.array(list(zip(_ord, np.ones_like(_ord, dtype=bool))), dtype=SeqOrder.INDEX_TYPE)

    def _update_positions(self):
        """
        An optimisation, whenever the positional state changes, this method must be called to
        maintain the current state in a separate array. This avoids unnecessary recalculation
        overhead.
        """
        # Masked sequences last, then by current position.
        sorted_indices = np.lexsort([self.order['pos'], ~self.order['mask']])
        for n, i in enumerate(sorted_indices):
            self.order[i]['pos'] = n
        self._positions = np.argsort(self.order['pos'])

    def remap_gapless(self, gapless_indices):
        """
        Recover the original, potentially sparse (gapped) indices from a dense (gapless) set
        of indices. Gaps originate from sequences being masked in the order. External tools
        often expect and return dense indices. When submitting changes to the current order
        state, it is important to first apply this method and reintroduce any gaps.

        Both a list/array of indices or a INDEX_TYPE array can be passed.

        :param gapless_indices: dense list of indices or an ndarray of type INDEX_TYPE
        :return: remappped indices with gaps (of a similar type to input)
        """
        # not as yet verified but this method is being replaced by the 50x faster numpy
        # alternative below. The slowless shows for large problems and repeated calls.
        # we ~could~ go further and maintain the shift array but this will require
        # consistent and respectful (fragile) use of mutator methods and not direct access on mask
        # or an observer.

        # the accumulated shifts due to masked sequences (the gaps).
        # we remove the masked sequences to make this array gapless
        shift = np.cumsum(~self.order['mask'])[self.order['mask']]

        # now reintroduce the gaps to the gapless representation supplied

        # handle our local type
        if isinstance(gapless_indices, np.ndarray) and gapless_indices.dtype == SeqOrder.INDEX_TYPE:
            remapped = []
            for oi in gapless_indices:
                remapped.append((oi['index'] + shift[oi['index']], oi['ori']))
            return np.array(remapped, dtype=SeqOrder.INDEX_TYPE)

        # handle plain collection
        else:
            remapped = []
            for oi in gapless_indices:
                remapped.append(oi + shift[oi])
            return np.array(remapped)

    def accepted_positions(self, copy=True):
        """
        The current positional order of only those sequences which have not been excluded by the mask.
        :param copy: return a copy
        :return: all accepted positons in order of index
        """
        return self.all_positions(copy=copy)[:self.count_accepted()]

    def all_positions(self, copy=True):
        """
        The current positional order of all sequences. Internal logic relegates masked sequences to always come
        last and ascending surrogate id order.

        :param copy: return a copy of the positions
        :return: all positions in order of index, masked or not.
        """
        if copy:
            _p = self._positions.copy()
        else:
            _p = self._positions
        return _p


    def gapless_positions(self):
        """
        A dense index range representing the current positional order without masked sequences. Therefore
        the returned array does not contain surrogate ids, but rather the relative positions of unmasked
        sequences, when all masked sequences have been discarded.

        :return: a dense index range of positional order, once all masked sequences have been discarded.
        """
        # accumulated shift from gaps
        gap_shift = np.cumsum(~self.order['mask'])
        # just unmasked sequences
        _p = np.argsort(self.order['pos'])
        _p = _p[:self.count_accepted()]
        # removing gaps leads to a dense range of indices
        _p -= gap_shift[_p]
        return _p

    def set_mask_only(self, _mask):
        """
        Set the mask state of all sequences, where indices in the mask map to
        sequence surrogate ids.

        :param _mask: mask array or list, boolean or 0/1 valued
        """
        _mask = np.asarray(_mask, dtype=bool)
        assert len(_mask) == len(self.order), 'supplied mask must be the same length as existing order'
        assert np.all((_mask == SeqOrder.ACCEPTED) | (_mask == SeqOrder.EXCLUDED)), \
            'new mask must be {} or {}'.format(SeqOrder.ACCEPTED, SeqOrder.EXCLUDED)

        # assign mask
        self.order['mask'] = _mask
        self._update_positions()


    def set_order_and_orientation(self, _ord, implicit_excl=False):
        """
        Set only the order, while ignoring orientation. An ordering is defined
        as a 1D array of the structured type INDEX_TYPE, where elements are the
        position and orientation of each indices.

        NOTE: This definition can be the opposite of what is returned by some
        ordering methods, and np.argsort(_v) should inverse the relation.

        NOTE: If the order includes only active sequences, setting implicit_excl=True
        the method will implicitly assume unmentioned ids are those currently
        masked. An exception is raised if a masked sequence is included in the order.

        :param _ord: 1d ordering
        :param implicit_excl: implicitly extend the order to include unmentioned excluded sequences.
        """
        assert _ord.dtype == SeqOrder.INDEX_TYPE, 'Wrong type supplied, _ord should be of INDEX_TYPE'

        if len(_ord) < len(self.order):
            # some sanity checks
            assert implicit_excl, 'Use implicit_excl=True for automatic handling ' \
                                  'of orders only mentioning accepted sequences'
            assert len(_ord) == self.count_accepted(), 'new order must mention all ' \
                                                       'currently accepted sequences'
            # those surrogate ids mentioned in the order
            mentioned = set(_ord['index'])
            assert len(mentioned & set(self.excluded())) == 0, 'new order and excluded must not ' \
                                                               'overlap when using implicit assignment'
            assert len(mentioned ^ set(self.accepted())) == 0, 'incomplete new order supplied,' \
                                                               'missing accepted ids'
            # assign the new orders
            self.order['pos'][_ord['index']] = np.arange(len(_ord), dtype=int)
            self.order['ori'][_ord['index']] = _ord['ori']
            # mask remaining, unmentioned indices
            _mask = np.zeros_like(self.mask_vector(), dtype=bool)
            _mask[_ord['index']] = True
            self.set_mask_only(_mask)
        else:
            # just a simple complete order update
            assert len(_ord) == len(self.order), 'new order was a different length'
            assert len(set(_ord['index']) ^ set(self.accepted())) == 0, 'incomplete new order supplied,' \
                                                                        'missing accepted ids'
            self.order['pos'][_ord['index']] = np.arange(len(_ord), dtype=int)
            self.order['ori'][_ord['index']] = _ord['ori']

        self._update_positions()

    

    def mask_vector(self):
        """
        :return: the current mask vector
        """
        return self.order['mask']

    def mask(self, _id):
        """
        Mask an individual sequence by its surrogate id

        :param _id: the surrogate id of sequence
        """
        self.order[_id]['mask'] = False
        self._update_positions()

    def count_accepted(self):
        """
        :return: the current number of accepted (unmasked) sequences
        """
        return self.order['mask'].sum()


    def accepted(self):
        """
        :return: the list surrogate ids for currently accepted sequences
        """
        return np.where(self.order['mask'])[0]


    def lengths(self, exclude_masked=False):
        # type: (bool) -> np.ndarray
        """
        Sequence lengths

        :param exclude_masked: True include only umasked sequencces
        :return: the lengths of sequences
        """
        if exclude_masked:
            return self.order['length'][self.order['mask']]
        else:
            return self.order['length']


class ContactMatrix:

    def __init__(self, bam_file, enzymes, seq_file, path ,metacc_folder, bin3c_folder, min_insert, bin_size, min_mapq_metacc, min_len_metacc, min_match_metacc,min_signal_metacc,min_mapq_bin3c ,min_len_bin3c,min_match_bin3c,min_signal_bin3c):

        ########input the parameter################
        ########################################################################
        ########################################################################
        #bam_file: alignment info of Hi-C library on contigs in bam#
        #enzymes: name of restriction enzymes used in Hi-C experiments#
        #seq_file: store the assembly contigs in fasta#
        #path: output path#
        #min_mapq_metacc: minimal mapping quality(default 30)#
        #min_len_metacc: minimal length of contigs(default 1000bp)
        #min_match_metacc: minimal match in cigar string (default 30)
        #min_signal_metacc : minimal cross contig signal (default 2)
        
        
        #min_mapq_bin3c = default 60
        #min_len_bin3c = default 1000
        #min_match_bin3c = default 10
        #min_signal_bin3c = default 5
        
        
        self.bam_file = bam_file
        self.enzymes = enzymes
        self.seq_file = seq_file
        self.path = path
        self.metacc_folder = metacc_folder
        self.bin3c_folder = bin3c_folder
        
        
        '''
        metacc parameters
        '''
        
        self.min_mapq_metacc = min_mapq_metacc
        self.min_len_metacc = min_len_metacc
        self.min_match_metacc = min_match_metacc
        self.min_signal_metacc = min_signal_metacc
        
        '''
        bin3c parameters
        '''
        
        self.min_mapq_bin3c = min_mapq_bin3c
        self.min_len_bin3c = min_len_bin3c
        self.min_match_bin3c = min_match_bin3c
        self.min_signal_bin3c = min_signal_bin3c
        self.order = None
        self.grouping = None 
        
        #fasta_info store the info from fasta file#
        #seq_info store the information of contigs from bam file#
        #seq_map store the contact map#
        self.fasta_info = {}
        
        '''
        metacc and bin3c have their own seq_info and seq_map
        '''
        self.seq_info_metacc = []
        self.seq_info_bin3c = []
        self.seq_map_metacc = None
        self.seq_map_bin3c = None
        
        '''
        primary_acceptance_mask is only for bin3c
        '''
        
        self.primary_acceptance_mask = None
        
        self.processed_map = None
        self.total_reads = None
        #the following min_insert and bin_size are only used by bin3C
        self.min_insert = min_insert
        self.bin_size = bin_size

        
        logger.info('Reading fasta file...')
        with open_input(seq_file) as multi_fasta:
            # prepare the site counter for the given experimental conditions
            # fasta_info is a dictionary of preparation of seq_info and seq_info is the true results
            site_counter = SiteCounter(enzymes,  is_linear=True)
            # get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count , desc='Analyzing contigs in reference fasta'):
                
                if len(seqrec) < min(min_len_metacc, min_len_bin3c):
                    continue
                self.fasta_info[seqrec.id] = {'sites': site_counter.count_sites(seqrec.seq),
                                         'length': len(seqrec)}

        logger.debug('There are {} contigs in reference fasta'.format(fasta_count))

        
        # now parse the header information from bam file
        ###########input the bam data###############
        #########seq_info contain the global contig information and we don't change the seq_info after being created##############
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
        ##test that BAM file is the correct sort order
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            # determine the set of active sequences
            # where the first filtration step is by length
            logger.info('Filtering metacc contigs by minimal length({})...'.format(self.min_len_metacc))
            logger.info('Filtering bin3c contigs by minimal length({})...'.format(self.min_len_bin3c))
            ref_count = {'seq_missing_metacc': 0, 'too_short_metacc': 0,'seq_missing_bin3c': 0, 'too_short_bin3c': 0}
            offset_metacc = 0
            localid_metacc = 0
            offset_bin3c = 0
            localid_bin3c = 0
            
            for n, (rname, rlen) in enumerate(zip(bam.references, bam.lengths)):
            # minimum length threshold
              if rlen >= min_len_metacc:
                
                try:
                    fa = self.fasta_info[rname]
                    assert fa['length'] == rlen, 'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)
                    self.seq_info_metacc.append(SeqInfo_metacc(localid_metacc , n , rname, fa['sites'], rlen, 0))
                    localid_metacc = localid_metacc + 1
                    offset_metacc += rlen
                except KeyError:
                    logger.info('Contig "{}" was not present in reference fasta'.format(rname))
                    ref_count['seq_missing_metacc'] += 1
                    
                 ######initially set the covcc coverage to be zero
                
              else:
                
                ref_count['too_short_metacc'] += 1
                
                
              if rlen >= min_len_bin3c:
                try:
                    fa = self.fasta_info[rname]
                    assert fa['length'] == rlen, 'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)
                    
                    
                    self.seq_info_bin3c.append(SeqInfo_bin3c(offset_bin3c, n, rname, rlen, fa['sites']))
                    ######initially set the covcc coverage to be zero
                
                    
                    localid_bin3c = localid_bin3c + 1
                    offset_bin3c += rlen
                except KeyError:
                    logger.info('Contig "{}" was not present in reference fasta'.format(rname))
                    ref_count['seq_missing_bin3c'] += 1
                           
              else:
                
                ref_count['too_short_bin3c'] += 1
                
                
              

            ####### total length of contigs##########
            ####### total_seq is number of contigs###
            self.total_len_metacc = offset_metacc
            self.total_seq_metacc = localid_metacc
            
            self.total_len_bin3c = offset_bin3c
            self.total_seq_bin3c = localid_bin3c
            

            del self.fasta_info

            if self.total_seq_metacc == 0:
                raise ImportError('No sequences in BAM file can be found in FASTA file for MetaCC')
            
            if self.total_seq_bin3c == 0:
                raise ImportError('No sequences in BAM file can be found in FASTA file for bin3C')
            
            logger.debug('{} contigs miss and {} contigs are too short for metacc'.format(ref_count['seq_missing_metacc'] , ref_count['too_short_metacc']))
            logger.debug('{} contigs miss and {} contigs are too short for bin3c'.format(ref_count['seq_missing_bin3c'] , ref_count['too_short_bin3c']))
            
            logger.debug('Accepted {} contigs covering {} bp for metacc'.format(self.total_seq_metacc, self.total_len_metacc))
            logger.debug('Accepted {} contigs covering {} bp for bin3c'.format(self.total_seq_bin3c, self.total_len_bin3c))
            logger.info('Counting reads in bam file...')

 
            self.total_reads = bam.count(until_eof=True)
            logger.debug('BAM file contains {0} alignments'.format(self.total_reads))
            
            # initialise the order
            self.order = SeqOrder(self.seq_info_bin3c)

            logger.info('Handling the alignments...')
            self._bin_map(bam)
            
            '''
            The primary_acceptance_mask is used by bin3c to do filtering
            '''
            self.set_primary_acceptance_mask()
                
        
        
        
        
        '''  
        
        following block is for metacc
        '''
        logger.info('Filtering contigs according to minimal signal of metacc({})...'.format(self.min_signal_metacc))
        contig_id = self.metacc_max_offdiag()
        logger.debug('{} contigs remain for metacc'.format(len(contig_id)))
        
        self.seq_map_metacc = self.seq_map_metacc.tolil()
        seq_temp_metacc = [] ###temporately store the sequence information#######
        
        for i , idn in enumerate(contig_id):
            seq = self.seq_info_metacc[idn]
            assert seq.localid_metacc == idn, 'the local index does not match the contact matrix index for metacc'
            seq_temp_metacc.append(SeqInfo_metacc(i , seq.refid , seq.name , seq.sites , seq.length , self.seq_map_metacc[idn , idn]))

        self.seq_info_metacc = seq_temp_metacc
        del seq_temp_metacc

        self.seq_map_metacc = self.seq_map_metacc.tocsr()
        self.seq_map_metacc = self.seq_map_metacc[contig_id , :]
        self.seq_map_metacc = self.seq_map_metacc.tocsc()
        self.seq_map_metacc = self.seq_map_metacc[: , contig_id]
        self.seq_map_metacc = self.seq_map_metacc.tocoo()
        del contig_id
        
        assert self.seq_map_metacc.shape[0] == len(self.seq_info_metacc), 'Filter error for metacc'
        
        ########change the diaganol entries of Hi-C matrix to zero#######
        self.seq_map_metacc = self.seq_map_metacc.tolil()
        self.seq_map_metacc.setdiag(0)
        
        #######Compute the Hi-C signals for each contig########
        self.seq_map_metacc = self.seq_map_metacc.tocsr()
        self.row_sum = np.matrix.tolist(self.seq_map_metacc.sum(axis=0))[0]
        self._write_contig_info_metacc()
        
        
        # create an initial acceptance mask
        
            
            
            
        
    def _bin_map(self, bam):
        """
        Accumulate read-pair observations from the supplied BAM file.
        Maps are initialized here. Logical control is achieved through initialisation of the
        ContactMap instance, rather than supplying this function arguments.

        :param bam: this instance's open bam file.
        """
        import tqdm

        def _simple_match_metacc(r):
            return r.mapping_quality >= _mapq_metacc

        def _strong_match_metacc(r):
            if r.mapping_quality < _mapq_metacc or r.cigarstring is None:
                return False
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            return cig[0] == 0 and cig[1] >= self.min_match_metacc

        # set-up match call
        _matcher_metacc = _strong_match_metacc if self.min_match_metacc else _simple_match_metacc
        
        def _simple_match_bin3c(r):
            return r.mapping_quality >= _mapq_bin3c

        def _strong_match_bin3c(r):
            if r.mapping_quality < _mapq_bin3c or r.cigarstring is None:
                return False
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            return cig[0] == 0 and cig[1] >= self.min_match_bin3c

        # set-up match call
        _matcher_bin3c = _strong_match_bin3c if self.min_match_bin3c else _simple_match_bin3c

        def next_informative(_bam_iter, _pbar):
            while True:
                r = next(_bam_iter)
                _pbar.update()
                if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
                    return r
       
        _seq_map_metacc = Sparse2DAccumulator(self.total_seq_metacc)
        _seq_map_bin3c = Sparse2DAccumulator(self.total_seq_bin3c)
        

        with tqdm.tqdm(total=self.total_reads) as pbar:

            # locals for read filtering
            _mapq_metacc = self.min_mapq_metacc
            _mapq_bin3c = self.min_mapq_bin3c

            _idx_metacc, _idx_bin3c = self.make_reverse_index('refid') #from global index to local index#
            _len = bam.lengths
     
            counts_metacc = OrderedDict({
                'accepted map_different_contig pairs for metacc': 0,
                'accepted map_same_contig pairs for metacc': 0,
                'ref_excluded pairs for metacc': 0,
                'poor_match pairs for metacc': 0,
                'single read for metacc':0})
            
            counts_bin3c = OrderedDict({
                'accepted map_different_contig pairs for bin3c': 0,
                'accepted map_same_contig pairs for bin3c': 0,
                'ref_excluded pairs for bin3c': 0,
                'poor_match pairs for bin3c': 0,
                'single read for bin3c':0})

            bam.reset()
            bam_iter = bam.fetch(until_eof=True)
            self.index1 = 0
            while True:
                self.index1 += 1
                try:
                    r1 = next_informative(bam_iter, pbar)
                    while True:
                        # read records until we get a pair
                        r2 = next_informative(bam_iter, pbar)
                        if r1.query_name == r2.query_name:
                            break
                        r1 = r2 ###if we don't get a pair, next _bam_iter
                        counts_metacc['single read for metacc'] += 1
                        counts_bin3c['single read for bin3c'] += 1
                        
                except StopIteration:
                    break

                if r1.reference_id in _idx_metacc and r2.reference_id in _idx_metacc:
                    if _matcher_metacc(r1) and _matcher_metacc(r2):
                        if r1.reference_id == r2.reference_id:
                            counts_metacc['accepted map_same_contig pairs for metacc'] += 1
                        else:
                            counts_metacc['accepted map_different_contig pairs for metacc'] += 1
                            
                        ix1_metacc = _idx_metacc[r1.reference_id]
                        ix2_metacc = _idx_metacc[r2.reference_id]
    
                        # maintain just a half-matrix
                        if ix2_metacc < ix1_metacc:
                            ix1_metacc, ix2_metacc = ix2_metacc, ix1_metacc
    
                        ix_metacc = (ix1_metacc , ix2_metacc)
                        if _seq_map_metacc.getitem(ix_metacc):
                            temp_value = _seq_map_metacc.getitem(ix_metacc) + 1
                            _seq_map_metacc.setitem(ix_metacc , temp_value)
                        else:
                            _seq_map_metacc.setitem(ix_metacc , 1)
                    
                    else:
                        counts_metacc['poor_match pairs for metacc'] += 1
                else:
                    counts_metacc['ref_excluded pairs for metacc'] += 1
                
                
                
                if r1.reference_id in _idx_bin3c and r2.reference_id in _idx_bin3c:
                    if _matcher_bin3c(r1)and _matcher_bin3c(r2):
                        if r1.reference_id == r2.reference_id:
                            counts_bin3c['accepted map_same_contig pairs for bin3c'] += 1
                        else:
                            counts_bin3c['accepted map_different_contig pairs for bin3c'] += 1
                        ix1_bin3c = _idx_bin3c[r1.reference_id]
                        ix2_bin3c = _idx_bin3c[r2.reference_id]
    
                        # maintain just a half-matrix
                        if ix2_bin3c < ix1_bin3c:
                            ix1_bin3c, ix2_bin3c = ix2_bin3c, ix1_bin3c
    
                        ix_bin3c = (ix1_bin3c , ix2_bin3c)
                        if _seq_map_bin3c.getitem(ix_bin3c):
                            temp_value = _seq_map_bin3c.getitem(ix_bin3c) + 1
                            _seq_map_bin3c.setitem(ix_bin3c , temp_value)
                        else:
                            _seq_map_bin3c.setitem(ix_bin3c , 1)
                    else:
                        counts_bin3c['poor_match pairs for bin3c'] += 1
                    
                else:
                  counts_bin3c['ref_excluded pairs for bin3c'] += 1
                  

                
        self.seq_map_metacc = _seq_map_metacc.get_coo()
        self.seq_map_bin3c = _seq_map_bin3c.get_coo()
        
        del _seq_map_metacc, _seq_map_bin3c, r1, r2, _idx_metacc, _idx_bin3c

        logger.debug('Pair accounting for metacc: {}'.format(counts_metacc))
        logger.debug('Pair accounting for bin3c: {}'.format(counts_bin3c))


    def make_reverse_index(self, field_name):
        """
        Make a reverse look-up (dict) from the chosen field in seq_info_metacc and seq_Info_bin3c to the internal index value
        of the given sequence. Non-unique fields will raise an exception.

        :param field_name: the seq_info field to use as the reverse.
        :return: internal array index of the sequence
        """
        rev_idx_metacc = {}
        rev_idx_bin3c = {}
        for n, seq in enumerate(self.seq_info_metacc):
            fv = getattr(seq, field_name)
            if fv in rev_idx_metacc:
                raise RuntimeError('field contains non-unique entries, a 1-1 mapping cannot be made for metacc')
            rev_idx_metacc[fv] = n
        
        for n, seq in enumerate(self.seq_info_bin3c):
            fv = getattr(seq, field_name)
            if fv in rev_idx_bin3c:
                raise RuntimeError('field contains non-unique entries, a 1-1 mapping cannot be made for bin3c')
            rev_idx_bin3c[fv] = n
        return rev_idx_metacc, rev_idx_bin3c


    def _write_contig_info_metacc(self):
        # Write detailed contig information
        with open(os.path.join(self.metacc_folder, 'tmp', 'contig_info_metacc.csv'), 'w') as out_tmp:
            out_tmp.write("name,sites,length,covcc,signal\n")
            for i, seq in enumerate(self.seq_info_metacc):
                out_tmp.write(f"{seq.name},{seq.sites},{seq.length},{seq.covcc},{self.row_sum[i]}\n")

        with open(os.path.join(self.metacc_folder , 'contig_info_metacc.csv'),'w') as out:
            out.write("name,sites,length\n")
            for seq in self.seq_info_metacc:
                out.write(str(seq.name)+ ',' + str(seq.sites)+ ',' + str(seq.length))
                out.write('\n')
        
        
    def metacc_max_offdiag(self):
        """
        Determine the maximum off-diagonal values of a given symmetric matrix. As this
        is assumed to be symmetric, we consider only the rows.

        :param _m: a scipy.sparse matrix
        :return: the off-diagonal maximum values
        """
        _m = self.seq_map_metacc
        assert scisp.isspmatrix(_m), 'Input matrix is not a scipy.sparse object'
        _m = _m.tolil(True)
        _diag = _m.tocsr().diagonal()
        _m.setdiag(0)
        _sig = np.asarray(_m.tocsr().max(axis=0).todense()).ravel()
        _contig_id = []
        for i in range(_m.shape[0]):
            if _sig[i] >= self.min_signal_metacc and _diag[i]>0 and self.seq_info_metacc[i].sites>0:
                _contig_id.append(i)
        del _m
        return _contig_id
    
    def get_primary_acceptance_mask(self):
        assert self.primary_acceptance_mask is not None, 'Primary acceptance mask has not be initialized'
        return self.primary_acceptance_mask.copy()

    def set_primary_acceptance_mask(self, min_len=None, min_sig=None, max_fold=None, update=False):
        """
        Determine and set the filter mask using the specified constraints across the entire
        contact map. The mask is True when a sequence is considered acceptable wrt to the
        constraints. The mask is also returned by the function for convenience.

        :param min_len: override instance value for minimum sequence length
        :param min_sig: override instance value for minimum off-diagonal signal (counts)
        :param max_fold: maximum locally-measured fold-coverage to permit
        :param update: replace the current primary mask if it exists
        :return: an acceptance mask over the entire contact map
        """
        assert max_fold is None, 'Filtering on max_fold is currently disabled'

        # If parameter based critiera were unset, use instance member values set at instantiation time
        if not min_len:
            min_len = self.min_len_bin3c
        if not min_sig:
            min_sig = self.min_signal_bin3c

        assert min_len, 'Filtering criteria min_len is None for bin3c'
        assert min_sig, 'Filtering criteria min_sig is None for bin3c'

        logger.debug('Setting primary acceptance mask with '
                     'filtering criterion min_len: {} min_sig: {} for bin3c'.format(min_len, min_sig))

        # simply return the current mask if it has already been determined
        # and an update is not requested
        if not update and self.primary_acceptance_mask is not None:
            logger.debug('Using existing mask')
            return self.get_primary_acceptance_mask()

        acceptance_mask = np.ones(self.total_seq_bin3c, dtype=bool)

        # mask for sequences shorter than limit
        _mask = self.order.lengths() >= min_len
        logger.debug('Minimum length threshold removing for bin3c: {}'.format(self.total_seq_bin3c - _mask.sum()))
        acceptance_mask &= _mask

        # mask for sequences weaker than limit
        
        
        signal = sparse_utils.max_offdiag(self.seq_map_bin3c)
        _mask = signal >= min_sig
        logger.debug('Minimum signal threshold removing for bin3c: {}'.format(self.total_seq_bin3c - _mask.sum()))
        acceptance_mask &= _mask

        # retain the union of all masks.
        self.primary_acceptance_mask = acceptance_mask

        logger.debug('Accepted sequences for bin3c: {}'.format(self.primary_acceptance_mask.sum()))

        return self.get_primary_acceptance_mask()

    def prepare_seq_map(self, norm=True, bisto=True, mean_type='geometric'):
        """
        Prepare the sequence map (seq_map_bin3c) by application of various filters and normalisations.

        :param norm: normalisation by sequence lengths
        :param bisto: make the output matrix bistochastic
        :param mean_type: when performing normalisation, use "geometric, harmonic or arithmetic" mean.
        """

        logger.info('Preparing sequence map for bin3c with full dimensions: {}'.format(self.seq_map_bin3c.shape))

        _mask = self.get_primary_acceptance_mask()

        self.order.set_mask_only(_mask)

        if self.order.count_accepted() < 1:
            raise NoneAcceptedException()

        _map = self.seq_map_bin3c.astype(float)

        # apply length normalisation if requested
        if norm:
            _map = self._norm_seq(_map, mean_type=mean_type, use_sites=True)
            logger.debug('Map normalized')

        # make map bistochastic if requested
        if bisto:
            # TODO balancing may be better done after compression
            _map, scl = self._bisto_seq(_map)
            # retain the scale factors
            self.bisto_scale = scl
            logger.debug('Map balanced')

        # cache the results for optional quick access
        self.processed_map = _map

    def get_subspace(self, permute=False, external_mask=None, marginalise=False, flatten=True,
                     dtype=float):
        """
        Using an already normalized full seq_map, return a subspace as indicated by an external
        mask or if none is supplied, the full map without filtered elements.

        The supplied external mask must refer to all sequences in the map.

        :param permute: reorder the map with the current ordering state
        :param external_mask: an external mask to combine with the existing primary mask
        :param marginalise: Assuming 4D NxNx2x2 tensor, sum 2x2 elements to become a 2D NxN
        :param flatten: convert a NxNx2x2 tensor to a 2Nx2N matrix
        :param dtype: return map with specific element type
        :return: subspace map
        """
        assert (not marginalise and not flatten) or np.logical_xor(marginalise, flatten), \
            'marginalise and flatten are mutually exclusive'

        # starting with the normalized map
        _map = self.processed_map.astype(dtype)

        # from a union of the sequence filter and external mask
        if external_mask is not None:
            _mask = self.get_primary_acceptance_mask()
            logger.info('Beginning with sequences after primary filtering: {}'.format(_mask.sum()))
            _mask &= external_mask
            logger.info('Active sequences after applying external mask: {}'.format(_mask.sum()))
            self.order.set_mask_only(_mask)

        # remove masked sequences from the map
        if self.order.count_accepted() < self.total_seq_bin3c:
            
            _map = sparse_utils.compress(_map.tocoo(), self.order.mask_vector())
            logger.info('After removing filtered sequences map dimensions: {}'.format(_map.shape))


        if permute:
            _map = self._reorder_seq(_map, flatten=flatten)
            logger.debug('Map reordered')

        return _map
    
    def _reorder_seq(self, _map, flatten=False):
        """
        Reorder a simple sequence map using the supplied map.

        :param _map: the map to reorder
        
        :return: ordered map
        """
        assert sp.isspmatrix(_map), 'reordering expects a sparse matrix type'

        _order = self.order.gapless_positions()
        

        assert _map.shape[0] == _order.shape[0], 'supplied map and unmasked order are different sizes'
        p = sp.lil_matrix(_map.shape)
        for i in range(len(_order)):
            p[i, _order[i]] = 1.
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _bisto_seq(self, _map):
        """
        Make a contact map bistochastic. This is another form of normslisation. Automatically
        handles 2D and 4D maps.

        :param _map: a map to balance (make bistochastic)
        :return: the balanced map
        """
        logger.debug('Balancing contact map')

        
        _map, scl = sparse_utils.kr_biostochastic(_map)
        return _map, scl

    def _get_sites(self):
        _sites = np.array([si.sites for si in self.seq_info_bin3c], dtype=float)
        # all sequences are assumed to have a minimum of 1 site -- even if not observed
        # TODO test whether it would be more accurate to assume that all sequences are under counted by 1.
        _sites[np.where(_sites == 0)] = 1
        return _sites

    def _norm_seq(self, _map, use_sites=True, mean_type='geometric'):
        """
        Normalise a simple sequence map in place by the geometric mean of interacting contig pairs lengths.
        The map is assumed to be in starting order.

        :param _map: the target map to apply normalisation
        :param use_sites: normalise matrix counts using observed sites, otherwise normalise
        using sequence lengths as a proxy
        :param mean_type: for length normalisation, choice of mean (harmonic, geometric, arithmetic)
        :return: normalized map
        """
        if use_sites:
            logger.debug('Doing site based normalisation')
            _sites = self._get_sites()
            _map = _map.astype(float)
            fast_norm_fullseq_bysite(_map.row, _map.col, _map.data, _sites)

        else:
            logger.debug('Doing length based normalisation')
        
            _mean_func = mean_selector(mean_type)
            _len = self.order.lengths().astype(float)
            _map = _map.tolil().astype(float)
            for i in range(_map.shape[0]):
                _map[i, :] /= np.fromiter((1e-3 * _mean_func(_len[i],  _len[j])
                                           for j in range(_map.shape[0])), dtype=float)
            _map = _map.tocsr()

        return _map
    
    def get_fields():
        """
        :return: the list of fields used in seq_info dict.
        """
        return SeqInfo_bin3c._fields


    def map_weight(self):
        """
        :return: the total map weight (sum ij)
        """
        return self.seq_map_bin3c.sum()

    def is_empty(self):
        """
        :return: True if the map has zero weight
        """
        return self.map_weight() == 0
