'''
Script for bin refinement using metahit
'''
import os
import sys
import logging
import subprocess
from MetaCC.Script.utils import save_object, gen_bins

if __name__ == '__main__':
    # Get the directory where the current script is located
    script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)-8s | %(asctime)s | %(message)s'
    )
    logger = logging.getLogger()
    
    try:
        # Get input parameters - in a real implementation you'd parse these from command line
        # For now, these would be defined based on the output from 5a_binning.py
        metacc_folder = sys.argv[1]  # Path to metacc output folder
        bin3c_folder = sys.argv[2]   # Path to bin3c output folder
        imputecc_folder = sys.argv[3]  # Path to imputecc output folder
        metahit_folder = sys.argv[4]  # Path to output folder
        
        logger.info('Starting bin refinement with metahit')
        
        metahit = os.path.join(script_directory, 'bin_refinement.sh')
        os.system("chmod 777 " + metahit)
        
        metahitCmd = metahit + " -t 80 -o "+metahit_folder + " -A " +os.path.join(metacc_folder,"BIN")+" -B "+ os.path.join(imputecc_folder , 'FINAL_BIN') + " -C " + os.path.join(bin3c_folder,"fasta")
        logger.info("metahitCmd : " + metahitCmd)
        output = os.popen(metahitCmd).read()
        logger.info(output)
        
        logger.info('Bin refinement completed successfully')
        
        # Paths
        final_bin_dir = os.path.join(metahit_folder, 'metahit_50_10_bins')
        final_fasta = os.path.join(metahit_folder, 'metahit_50_10_bins.fa')
        cluster_txt = os.path.join(metahit_folder, 'cluster.txt')
        
        # Manually create the merged .fa file
        with open(final_fasta, 'w') as outfile:
            for fname in sorted(os.listdir(final_bin_dir)):
                if fname.endswith('.fa'):
                    with open(os.path.join(final_bin_dir, fname)) as infile:
                        outfile.write(infile.read())
        
        # Generate contig-to-bin mapping directly
        final_clustering = {}
        
        for fname in sorted(os.listdir(final_bin_dir)):
            if not fname.endswith(".fa"):
                continue
            bin_id = os.path.splitext(fname)[0]  # e.g., 'bin.1.fa' → 'bin.1'
            with open(os.path.join(final_bin_dir, fname)) as f:
                for line in f:
                    if line.startswith(">"):
                        contig = line[1:].strip()
                        final_clustering[contig] = bin_id
        
        # Save cluster.txt manually
        with open(cluster_txt, 'w') as f:
            for contig, bin_id in final_clustering.items():
                f.write(f"{contig}\t{bin_id}\n")
        
        # Save final_bins.p
        save_object(os.path.join(metahit_folder, 'final_bins.p'), final_clustering)
        
    except Exception as e:
        logger.error(f'Error during bin refinement: {str(e)}')
        sys.exit(1)