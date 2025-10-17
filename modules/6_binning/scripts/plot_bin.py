'''
libraries for metacc
'''

from MetaCC.Script.exceptions import ApplicationException
from MetaCC.Script.utils import load_object, save_object, make_dir, gen_bins, gen_sub_bins


#######Import python packages
import subprocess
import argparse
import warnings
import logging
import shutil
import sys
import os
import fileinput
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import scipy.sparse as scisp
import seaborn

##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '1.0.0, released at 10/2024'

def replace_in_file(file_path):
    with fileinput.input(file_path, inplace=True) as file:
        for line in file:
            if ">" and ":"  in line:
              new_line = line.split(' ')[1].split(':')[1]
              print(">"+new_line)
            else:
              print(line,end="")
              
              
def enable_clusters(contact_map, clustering, cl_list=None, ordered_only=True, min_extent=None):
    """
    Given a clustering and list of cluster ids (or none), enable (unmask) the related sequences in
    the contact map. If a requested cluster has not been ordered, it will be dropped.

    :param contact_map: an instance of ContactMap to cluster
    :param clustering: a clustering solution to the contact map
    :param cl_list: a list of cluster ids to enable or None (all ordered clusters)
    :param ordered_only: include only clusters which have been ordered
    :param min_extent: include only clusters whose total extent is greater
    :return: the filtered list of cluster ids in ascending numerical order
    """

    # start with all clusters if unspecified
    if cl_list is None:
        cl_list = list(clustering.keys())

    # use instance criterion if not explicitly set
    if min_extent is None:
        min_extent = contact_map.min_extent

    if min_extent:
        cl_list = [k for k in cl_list if clustering[k]['extent'] >= min_extent]
        logger.info('Clusters passing minimum extent criterion: {}'.format(len(cl_list)))
        if len(cl_list) == 0:
            raise NoRemainingClustersException(
                'No clusters passed min_extent criterion of >= {}'.format(min_extent))

    if ordered_only:
        # drop any clusters that have not been ordered
        cl_list = [k for k in cl_list if 'order' in clustering[k]]
        logger.info('Clusters passing ordered-only criterion: {}'.format(len(cl_list)))
        if len(cl_list) == 0:
            raise NoRemainingClustersException(
                'No clusters passed ordered-only criterion')

    # impose a consistent order
    cl_list = sorted(cl_list)

    if ordered_only:
        # use the determined order and orientation
        cmb_ord = np.hstack([clustering[k]['order'] for k in cl_list])
    else:
        # arbitrary order and orientation
        cmb_ord = np.hstack([SeqOrder.asindex(clustering[k]['seq_ids']) for k in cl_list])

    if len(cmb_ord) == 0:
        raise NoRemainingClustersException('No requested cluster contained ordering information')

    logger.info('Total number of sequences in the clustering: {}'.format(len(cmb_ord)))

    # prepare the mask
    _mask = np.zeros_like(contact_map.order.mask_vector(), dtype=np.bool)
    _mask[cmb_ord['index']] = True
    _mask &= contact_map.get_primary_acceptance_mask()
    logger.info('After joining with active sequence mask map: {}'.format(_mask.sum()))
    contact_map.order.set_mask_only(_mask)
    contact_map.order.set_order_and_orientation(cmb_ord, implicit_excl=True)

    return cl_list


def plot_clusters(contact_map, fname, clustering, cl_list=None, simple=True, permute=False, max_image_size=None,
                  ordered_only=False, min_extent=None, use_taxo=False, flatten=False, **kwargs):
    """
    Plot the contact map, annotating the map with cluster names and boundaries.

    For large contact maps, block reduction can be employed to reduce the size for plotting purposes. Using
    block_reduction=2 will reduce the map dimensions by a factor of 2. Must be integer.

    :param contact_map: an instance of ContactMap to cluster
    :param fname: output file name
    :param clustering: the cluster solution
    :param cl_list: the list of cluster ids to include in plot. If none, include all ordered clusters
    :param simple: True plot seq map, False plot the extent map
    :param permute: permute the map with the present order
    :param max_image_size:  maximum allowable image size before rescale occurs
    :param ordered_only: include only clusters which have been ordered
    :param min_extent: include only clusters whose total extent is greater
    :param use_taxo: use taxonomic information within clustering, assuming it exists
    :param flatten: for tip-based, flatten matrix rather than marginalise
    :param kwargs: additional options passed to plot()
    """

    if cl_list is None:
        logger.info('Plotting heatmap of complete solution')
    else:
        logger.info('Plotting heatmap for {} specified clusters'.format(len(cl_list)))


    # now build the list of relevant clusters and setup the associated mask
    cl_list = enable_clusters(contact_map, clustering, cl_list=cl_list, ordered_only=ordered_only,
                              min_extent=min_extent)

    
    tick_locs = np.cumsum([0] + [len(clustering[k]['seq_ids']) for k in cl_list])

    plot(contact_map,fname, permute=permute, simple=simple, tick_locs=tick_locs, 
                     max_image_size=max_image_size, flatten=flatten, **kwargs)
    

def plot(cm, fname, simple=False, tick_locs=None, tick_labs=None, norm=True, permute=False, pattern_only=False,
         dpi=180, width=25, height=22, zero_diag=None, alpha=0.01, robust=False, max_image_size=None,
         flatten=False):
    """
    Plot the contact map. This can either be as a sparse pattern (requiring much less memory but without visual
    cues about intensity), simple sequence or full binned map and normalized or permuted.

    :param fname: output file name
    :param tick_locs: major tick locations (minors take the midpoints)
    :param tick_labs: minor tick labels
    :param simple: if true, sequence only map plotted
    :param norm: normalize intensities by geometric mean of lengths
    :param permute: reorder map to current order
    :param pattern_only: plot only a sparse pattern (much lower memory requirements)
    :param dpi: adjust DPI of output
    :param width: plot width in inches
    :param height: plot height in inches
    :param zero_diag: set bright cm-interactions to zero
    :param alpha: log intensities are log (x + alpha)
    :param robust: use seaborn robust dynamic range feature
    :param max_image_size: maximum allowable image size before rescale occurs
    :param flatten: for tip-based, flatten matrix rather than marginalise
    """

    plt.style.use('ggplot')

    fig = plt.figure()
    fig.set_figwidth(width)
    fig.set_figheight(height)
    ax = fig.add_subplot(111)


    # prepare the map if not already done. This overwrites
    # any current ordering mask beyond the primary acceptance mask
    if cm.processed_map is None:
        cm.prepare_seq_map(norm=norm, bisto=True)
    _map = cm.get_subspace(permute=permute, marginalise=False if flatten else True, flatten=flatten)
    # unless requested, zero diagonal for simple plots as its intensity tends to obscure detail
    if zero_diag is None:
        _map.setdiag(0)
    # amplify values for plotting
    _map *= 100
   

    if pattern_only:
        # sparse matrix plot, does not support pixel intensity
        if zero_diag:
            _map.setdiag(0)
        ax.spy(_map.tocsr(), markersize=5 if simple else 1)

    else:
        # a dense array plot

        # if too large, reduced it while sparse.
        if max_image_size is not None:
            full_size = _map.shape
            if np.max(full_size) > max_image_size:
                reduce_factor = int(np.ceil(np.max(full_size) / float(max_image_size)))
                logger.info('Full {} image reduction factor: {}'.format(full_size, reduce_factor))
                # downsample the map
                _map = sparse_utils.downsample(_map, reduce_factor)
                # ticks adjusted to match
                tick_locs = np.floor(tick_locs.astype(np.float) / reduce_factor)
                logger.info('Map has been reduced from {} to {}'.format(full_size, _map.shape))

        _map = _map.toarray()

        if zero_diag:
            logger.debug('Removing diagonal')
            np.fill_diagonal(_map, 0)

        _map = np.log(_map + alpha)

        logger.debug('Making raster image')
        seaborn.heatmap(_map, robust=robust, square=True, linewidths=0, ax=ax, cbar=False)

    if tick_locs is not None:

        plt.tick_params(axis='both', which='both',
                        right=False, left=False, bottom=False, top=False,
                        labelright=False, labelleft=False, labelbottom=False, labeltop=False)

        if tick_labs is not None:
            min_labels = ticker.FixedFormatter(tick_labs)
            ax.tick_params(axis='y', which='minor', left=True, labelleft=True, labelsize=10)

            min_ticks = ticker.FixedLocator(tick_locs[:-1] + 0.5 * np.diff(tick_locs))

            ax.yaxis.set_minor_formatter(min_labels)
            ax.yaxis.set_minor_locator(min_ticks)

        # seaborn will not display the grid, so we make our own.
        ax.hlines(tick_locs, *ax.get_xlim(), color='grey', linewidth=0.5, linestyle='-.')
        ax.vlines(tick_locs, *ax.get_ylim(), color='grey', linewidth=0.5, linestyle='-.')

    logger.debug('Saving plot')
    fig.tight_layout()
    plt.savefig(fname, dpi=dpi)
    plt.close(fig)

if __name__ == '__main__':
    
    
    
    def mk_version():
        return 'Metahit v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg
    
    
    
    script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    
    
    global_parser.add_argument('--contact-map', help='contact map object')
    global_parser.add_argument('--BIN', help ='contact map object')
    global_parser.add_argument('--OUTDIR', help='Output directory')
    

    args = global_parser.parse_args()


    
    
    
    try:
        make_dir(args.OUTDIR)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
    
 
           
    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    log_path = os.path.join(args.OUTDIR, 'plot.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))
    
    
    try:
        
        logger.info('Loading contact maps from {}'.format(args.contact_map))
        cm = load_object(args.contact_map)
        
        #########Scan the marker gene to determine the hyperparameter in the Leiden clustering#########
        

        # cluster the entire map
        #clustering = cluster_map(cm, method='infomap', seed=args.seed, work_dir=bin3c_folder)
        # generate report per cluster

        # serialize full clustering object
        
        clustering = args.BIN
        

        plot_clusters(cm, os.path.join(args.OUTDIR, 'cluster_plot.png'), clustering,
                          max_image_size=args.max_image, ordered_only=False, simple=False, permute=True)


    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
