'''
libraries for metacc
'''
from raw_contact_both import ContactMatrix
from MetaCC.Script.normalized_contact import NormCCMap, NormCCMap_LC
from MetaCC.Script.predict_species_number import gen_bestk
from MetaCC.Script.cluster import ClusterBin
from MetaCC.Script.exceptions import ApplicationException
from MetaCC.Script.utils import load_object, save_object, make_dir, gen_bins, gen_sub_bins, make_random_seed
from MetaCC.Script.normcc import normcc
import scipy.sparse as scisp
import argparse
import warnings
import logging
import shutil
import sys
import os
import pandas as pd

'''
libraries for bin3c
'''

#!/usr/bin/env python
from bin3C_python3.mzd.cluster import *
from bin3C_python3.mzd.exceptions import ApplicationException
#from bin3C_python3.mzd.io_utils import load_object, save_object
from bin3C_python3.mzd.utils import *
import logging
import sys


'''
libraries for imputecc
'''

#########The structure of the main script is modified from bin3C########
from ImputeCC.imputation import ImputeMatrix
from ImputeCC.exceptions import ApplicationException
from ImputeCC.pre_common import get_bin_dirs
from ImputeCC.pre_profile import Profile
from ImputeCC.utility import save_object, load_object, make_dir, gen_bins
from ImputeCC.final_cluster import FinalCluster
from ImputeCC.pre_clustering import Clust4CheckM, PreCluster

#######Import python packages
import subprocess
import argparse
import warnings
import logging
import shutil
import sys
import os
import fileinput

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
    
    
    global_parser.add_argument('--FASTA', help='Reference fasta sequence')
    global_parser.add_argument('--BAM', help='Input bam file in query order')
    global_parser.add_argument('--OUTDIR', help='Output directory')
    global_parser.add_argument(
        '-e', '--enzyme', metavar='NEB_NAME', action='append', required=False, default=["Sau3AI,MluCI"],
        help='Comma-separated enzyme names (e.g., Sau3AI,MluCI). Default: Sau3AI,MluCI'
    )
    
    
    global_parser.add_argument('--metacc-min-len', type=int,
                           help='Minimum acceptable contig length [1000]')
    global_parser.add_argument('--metacc-min-signal', type=int,
                           help='Minimum acceptable Hi-C signal [2]')
    global_parser.add_argument('--metacc-min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    global_parser.add_argument('--metacc-min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    global_parser.add_argument('--thres', type=float,
                           help='the fraction of discarded NormCC-normalized Hi-C contacts [0.05]')
    global_parser.add_argument('--num-gene', type=int,
                               help='Number of maker genes detected, automatically detected if not input')
    
    global_parser.add_argument('--seed', type=int, default=None,
                               help='Random seed')
    global_parser.add_argument('--metacc-min-binsize', type=int,
                               help='Minimum bin size used in output [150000]')
    
    
    
    
    global_parser.add_argument('--bin3c-min-len', type=int,
                           help='Minimum acceptable contig length [1000]')
    global_parser.add_argument('--bin3c-min-signal', type=int,
                           help='Minimum acceptable Hi-C signal [2]')
    global_parser.add_argument('--bin3c-min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    global_parser.add_argument('--bin3c-min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    global_parser.add_argument('--bin3c-min-extent', type=int,
                               help='Minimum cluster extent used in output [50000]')
    
    
    global_parser.add_argument('--no-fasta', default=False, action='store_true',
                             help='Do not generate cluster FASTA files')
    
    
    global_parser.add_argument('--no-report', default=False, action='store_true',
                             help='Do not generate a cluster report')
    
    global_parser.add_argument('--no-spades', default=False, action='store_true',
                             help='Assembly was not done using SPAdes')
    
    global_parser.add_argument('--only-large', default=False, action='store_true',
                             help='Only write FASTA for clusters longer than min_extent')
    
    
    
    
    global_parser.add_argument('-t', '--threads', default=30, type=int, help="the number of threads. default is 30.")
    
    
    
    global_parser.add_argument('--gene-cov', type=float, help='gene coverage used in detecting marker genes, default 0.9')
    global_parser.add_argument('--rwr-rp', type=float, help='random walk restart probability, default 0.5')
    global_parser.add_argument('--rwr-thres', type=int, help='cut-off to maintain sparsity in each random walk step, default 80')
    global_parser.add_argument('--max-markers', type=int, help='maximum number of contigs with marker genes, default 8000')
    
    global_parser.add_argument('--intra', type=int, help='percentile threshold to assign the contigs to preliminary bins in pre-clustering step, default 50')
    global_parser.add_argument('--inter', type=int, help='percentile threshold to assign the contigs to new bins in pre-clustering step, default 0')
    
    global_parser.add_argument('--cont-weight', type=float, help='coefficient of completeness - cont_weight * completeness, default 2')
    global_parser.add_argument('--min-comp', type=float, help='minimum completeness of bin to consider during bin selection process, default 50')
    global_parser.add_argument('--max-cont' , type=float, help='maximum contamination of bin to consider during bin selection process, default 10')
    global_parser.add_argument('--report-quality', type=float, help="minimum quality of bin to report, default 10")
    global_parser.add_argument('--imputecc-min-binsize', type=int, help='Minimum bin size used in output [100000]')





    args = global_parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)
    
    
    
    
    # Create folders for metacc, bin3c and imputecc folder
    metacc_folder = os.path.join(args.OUTDIR , 'metacc')
    bin3c_folder = os.path.join(args.OUTDIR , 'bin3c')
    imputecc_folder = os.path.join(args.OUTDIR , 'imputecc')
    metacc_temp_folder = os.path.join(metacc_folder , 'tmp')
    metahit_folder = os.path.join(args.OUTDIR , 'metahit_bin_refinement')
    
    
    
    try:
        make_dir(args.OUTDIR)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
    
    try:
        make_dir(metacc_folder)
    except IOError:
        print('Error: cannot find out directory or the metacc directory already exists')
        sys.exit(1)
        
    try:
        make_dir(bin3c_folder)
    except IOError:
        print('Error: cannot find out directory or the bin3c directory already exists')
        sys.exit(1)
    
    try:
        make_dir(imputecc_folder)
    except IOError:
        print('Error: cannot find out directory or the imputecc directory already exists')
        sys.exit(1)
        
    
    
    
    if not os.path.exists(metacc_temp_folder):
        os.mkdir(metacc_temp_folder)
    else:
        #shutil.rmtree(metacc_temp_folder)           
        os.mkdir(metacc_temp_folder)
    
    
    imputecc_temp_folder = os.path.join(imputecc_folder , 'tmp')
    if not os.path.exists(imputecc_temp_folder):
        os.mkdir(imputecc_temp_folder)
    else:
        #shutil.rmtree(imputecc_temp_folder)           
        os.mkdir(imputecc_temp_folder)
    
    
    
    # Create Intermediate folder
    # Intermediate folder is not deleted
    interm_folder = os.path.join(imputecc_folder, 'intermediate')
    if not os.path.exists(interm_folder):
        os.mkdir(interm_folder)
        
    
        
    
           
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
    log_path = os.path.join(args.OUTDIR, 'bin_refinement.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))
    
    
    
    
    bin3c_runtime_defaults = {
        'min_reflen': 1000,
        'min_signal': 5,
        'min_extent': 50000,
        'min_mapq': 60,
        'strong': 10
    }
    
    metacc_runtime_defaults = {
        'min_len': 1000,
        'min_signal': 2,
        'min_mapq': 30,
        'min_match': 30,
        'thres': 0.05,
        'min_binsize':150000
    }
    
    imputecc_runtime_defaults = {
        'gene_cov': 0.9,
        'rwr_rp': 0.5,
        'rwr_thres': 80,
        'max_markers': 8000,
        'intra': 50,
        'inter': 0,
        'min_binsize':100000,
        'cont_weight': 2,
        'min_comp': 50.0,
        'max_cont': 10.0,
        'report_quality': 10.0
    }
    
    enzyme = []
    for e in args.enzyme:
        enzyme.extend(e.split(","))

    try:
        logger.info('Begin constructing raw contact matrix for metacc and bin3c...')
        cm = ContactMatrix(
            args.BAM,
            enzyme,
            args.FASTA,
            args.OUTDIR,
            metacc_folder,
            bin3c_folder,
            min_insert=0,
            bin_size=None,
            min_mapq_metacc=ifelse(args.metacc_min_mapq, metacc_runtime_defaults['min_mapq']),
            min_len_metacc=ifelse(args.metacc_min_len, metacc_runtime_defaults['min_len']),
            min_match_metacc=ifelse(args.metacc_min_match, metacc_runtime_defaults['min_match']),
            min_signal_metacc=ifelse(args.metacc_min_signal, metacc_runtime_defaults['min_signal']),
            min_mapq_bin3c=ifelse(args.bin3c_min_mapq, bin3c_runtime_defaults['min_mapq']),
            min_len_bin3c=ifelse(args.bin3c_min_len, bin3c_runtime_defaults['min_reflen']),
            min_match_bin3c=ifelse(args.bin3c_min_match, bin3c_runtime_defaults['strong']),
            min_signal_bin3c=ifelse(args.bin3c_min_signal, bin3c_runtime_defaults['min_signal'])
        )
        logger.info('Raw contact matrix construction finished for metacc and bin3c.')

        save_object(os.path.join(metacc_folder, 'contact_map.p'), cm)

        '''


          #The following is the normalization process of MetaCC

        '''

        logger.info('Begin normalizing raw contacts by NormCC for metacc raw matrix')
        
        contig_file = os.path.join(metacc_temp_folder , 'contig_info_metacc.csv')
        norm_result = normcc(contig_file)
        contig_info_df = pd.read_csv(contig_file)
        
        ######Construct normalized matrix of Hi-C interaction maps#############
        hzmap = NormCCMap(metacc_folder,
                        contig_info_df,
                        cm.seq_map_metacc,
                        norm_result,
                        thres = ifelse(args.thres, metacc_runtime_defaults['thres']))
                        
        logger.info('NormCC normalization finished')
                    
    
        #shutil.rmtree(metacc_temp_folder) ######Remove all intermediate files#######
        scisp.save_npz(os.path.join(metacc_folder, 'Normalized_contact_matrix.npz'), hzmap.seq_map.tocsr())
        save_object(os.path.join(metacc_folder, 'NormCC_normalized_contact'), hzmap)
        
        logger.info('MetaCC Normalization results have been saved')
        
        
        if not os.path.exists(os.path.join(metacc_folder , 'NormCC_normalized_contact.gz')):
            raise IOError('Please run the NormCC normalization step before binning')
        
        ###########Load the normalization instance to get access to the normalized Hi-C contact maps##########
        logger.info('Loading normalized contact maps by NormCC from: {}'.format(os.path.join(metacc_folder , 'NormCC_normalized_contact.gz')))
        hzmap = load_object(os.path.join(metacc_folder , 'NormCC_normalized_contact.gz'))
        
        #########Scan the marker gene to determine the hyperparameter in the Leiden clustering#########
        if args.num_gene is None:
            logger.info('Begin scanning marker genes...')
            args.num_gene = gen_bestk(metacc_folder , args.FASTA)
            if args.num_gene == 0:
                logger.warning('No marker gene is detected from the assembled contigs!')
            
        logger.info('There are {} marker genes in the assembled contigs'.format(args.num_gene))
        
        if not args.seed:
            args.seed = make_random_seed()
            
        logger.info('The random seed for clustering is {}'.format(args.seed))
            
        cluster_process = ClusterBin(metacc_folder , hzmap.name , hzmap.len , hzmap.seq_map ,
                                        ifelse(args.metacc_min_binsize, metacc_runtime_defaults['min_binsize']), args.num_gene, args.seed)
        logger.info('Writing bins...')
        gen_bins(args.FASTA , os.path.join(metacc_temp_folder , 'cluster.txt') , os.path.join(metacc_folder ,'BIN'))
        #shutil.rmtree(metacc_temp_folder) ######Remove all intermediate files#######
        logger.info('MetaCC binning fininshed.')


        
    
        '''
        #The following is the normalization and clustering process of bin3C
    
        '''
        
        
        cm.min_extent = ifelse(args.bin3c_min_extent, bin3c_runtime_defaults['min_extent'])

        
        

        # cluster the entire map
        clustering = cluster_map(cm, method='infomap', seed=args.seed, work_dir=bin3c_folder)
        # generate report per cluster
        cluster_report(cm, clustering, is_spades=not args.no_spades)
        
        # serialize full clustering object
        save_object(os.path.join(bin3c_folder, 'clustering.p'), clustering)
        
        
        if not args.no_report:
            # write a tabular report
            write_report(os.path.join(bin3c_folder, 'cluster_report.csv'), clustering)

        if not args.no_fasta:
            # write per-cluster fasta files, also separate ordered fasta if an ordering exists
            write_fasta(cm, bin3c_folder, clustering, source_fasta=args.FASTA, clobber=True,
                        only_large=args.only_large)



        

        '''
        #the following is for imputecc
        '''


        CONTIG_INFO = os.path.join(metacc_folder,"contig_info_metacc.csv")
        HIC_MATRIX = os.path.join(metacc_folder,"Normalized_contact_matrix.npz")


        logger.info('The ImputeCC binning pipeline starts...')
        
        marker_file = os.path.join(metacc_temp_folder , 'contigs.hmmout')
        if os.path.exists(marker_file):
            logger.info('Existing single-copy marker gene file is detected.')
        else:
            if args.num_gene is None:
                logger.info('Begin detecting single-copy marker genes from assembled contigs...')
                args.num_gene = gen_bestk(interm_folder , args.FASTA)
                if args.num_gene == 0:
                    logger.warning('No marker gene is detected from the assembled contigs!')      
                marker_file = os.path.join(interm_folder , 'contigs.hmmout')
                logger.info('There are {} marker genes in the assembled contigs'.format(args.num_gene))
            
              
        if not os.path.exists(marker_file):
            raise ModuleNotFoundError('Single-copy marker gene detection failed!')
        else:
            logger.info('Single-copy marker gene detection finished.')
        
        logger.info('Begin CRWR Imputation...')
        imp = ImputeMatrix(CONTIG_INFO,
                            HIC_MATRIX,
                            marker_file,
                            gene_cov = ifelse(args.gene_cov, imputecc_runtime_defaults['gene_cov']),
                            rwr_rp = ifelse(args.rwr_rp, imputecc_runtime_defaults['rwr_rp']),
                            rwr_thres= ifelse(args.rwr_thres, imputecc_runtime_defaults['rwr_thres']),
                            max_markers= ifelse(args.max_markers, imputecc_runtime_defaults['max_markers'])                                                  
                            )
        
        save_object(os.path.join(interm_folder, 'ImputeCC_storage'), imp)
        logger.info('CRWR Imputation finished.')
        
        #######Construct preliminary bins#####
        logger.info('Begin preclustering marker-gene-containing contigs...')
        bins, bin_of_contigs = PreCluster(imp.marker_contig_counts, imp.marker_contigs, imp.contig_markers,
                                         imp.imputed_matrix, imp.dict_contigRevLocal,
                                         intra= ifelse(args.intra , imputecc_runtime_defaults['intra']),
                                         inter = ifelse(args.inter , imputecc_runtime_defaults['inter']))
        logger.info('Preclustering finished with {} preliminary bins established.'.format(len(bins)))
        
        
        logger.info('Begin final binning for all contigs utilizing the information of preliminary bins...')            
        checkm_bac_gene_table = os.path.join(interm_folder, 'checkm_gene_table' , 'checkm_bac_gene_table.tsv')
        checkm_ar_gene_table = os.path.join(interm_folder, 'checkm_gene_table' , 'checkm_ar_gene_table.tsv')
        if os.path.exists(checkm_ar_gene_table) and os.path.exists(checkm_bac_gene_table):
            logger.info('CheckM lineage-specific gene tables detected.')
        else:
            logger.info('Detect lineage-specific genes...')
            Clust4CheckM(args.FASTA, CONTIG_INFO, HIC_MATRIX, imputecc_folder)
            output_dir = os.path.join(imputecc_temp_folder, 'out_checkm_qa')
            bin_file = os.path.join(imputecc_temp_folder , 'dir4checkm.tsv')
            cpus = args.threads
            bin_dirs = get_bin_dirs(bin_file)
            profile = Profile(cpus)
            profile.run(bin_dirs,
                        output_dir) 
            os.mkdir(os.path.join(interm_folder, 'checkm_gene_table'))

            mv1Cmd = 'mv ' + os.path.join(imputecc_temp_folder, 'out_checkm_qa', 'binning_methods', 'INITIAL_BIN', 
                                        'checkm_bac', 'marker_gene_table.tsv') + ' ' + checkm_bac_gene_table
            
            mv2Cmd = 'mv ' + os.path.join(imputecc_temp_folder, 'out_checkm_qa', 'binning_methods', 'INITIAL_BIN', 
                                        'checkm_ar', 'marker_gene_table.tsv') + ' ' + checkm_ar_gene_table
            os.system(mv1Cmd)
            os.system(mv2Cmd)
            logger.info('Lineage-specific gene detection finished.')
            
        cluster_process = FinalCluster(imp.contig_info, imp.contig_local, 
             imp.dict_contigRev, imp.dict_contigRevLocal, imp.dict_contig_len, 
             imp.contig_markers, bins, bin_of_contigs,
             imp.normcc_matrix, imp.imputed_matrix, 
             checkm_bac_gene_table, checkm_ar_gene_table, imputecc_folder,
             intra= ifelse(args.intra , imputecc_runtime_defaults['intra']),
             inter = ifelse(args.inter , imputecc_runtime_defaults['inter']),
             cont_weight = ifelse(args.cont_weight , imputecc_runtime_defaults['cont_weight']),
             min_comp = ifelse(args.min_comp , imputecc_runtime_defaults['min_comp']), 
             max_cont = ifelse(args.max_cont , imputecc_runtime_defaults['max_cont']), 
             report_quality = ifelse(args.report_quality , imputecc_runtime_defaults['report_quality']),
             min_binsize = ifelse(args.imputecc_min_binsize , imputecc_runtime_defaults['min_binsize']))
        
        logger.info('Final binning finished!')
        logger.info('Writing the final bins...')
        gen_bins(args.FASTA , os.path.join(imputecc_temp_folder , 'cluster_imputecc.txt') , os.path.join(imputecc_folder ,'FINAL_BIN'))
        #shutil.rmtree(imputecc_temp_folder) ######Remove all intermediate files#######
        logger.info('The ImputeCC binning pipeline fininshed.')
        
        
        
        logger.info('Change the naming of contigs for bin3C binning results') 
        
        for file in os.listdir(os.path.join(bin3c_folder,"fasta")):
            filename = file
            #print(filename)
            if filename.endswith(".fna"):
                 replace_in_file(os.path.join(os.path.join(bin3c_folder,"fasta"), filename))
         
        logger.info('Change the naming of contigs for bin3C binning results fininshed.')

    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
