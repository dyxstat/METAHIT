'''
Script for bin refinement using metahict
'''
import os
import sys
import logging
import subprocess
import argparse
import shutil
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
        parser = argparse.ArgumentParser(description="Integrate MetaCC, bin3C, and ImputeCC bins with METAHICT.")
        parser.add_argument("metacc_folder", help="MetaCC output folder")
        parser.add_argument("bin3c_folder", help="bin3C output folder")
        parser.add_argument("imputecc_folder", help="ImputeCC output folder")
        parser.add_argument("metahict_folder", help="METAHICT integration output folder")
        parser.add_argument("-t", "--threads", type=int, default=16, help="Number of CPU threads (default=16)")
        parser.add_argument("--min-completeness", type=float, default=50, help="Final METAHICT minimum completeness (default=50)")
        parser.add_argument("--max-contamination", type=float, default=10, help="Final METAHICT maximum contamination (default=10)")
        parser.add_argument("--contamination-penalty", type=float, default=5, help="Penalty used in completeness - penalty * contamination (default=5)")
        parser.add_argument("--min-input-bin-size", type=int, default=50000, help="Minimum input bin FASTA file size before refinement (default=50000 bytes)")
        parser.add_argument("--max-input-bin-size", type=int, default=20000000, help="Maximum input bin FASTA file size before refinement (default=20000000 bytes)")
        parser.add_argument("--binning-refiner-min-size", type=int, default=524288, help="Minimum refined bin size for Binning_refiner (default=524288 bp)")
        parser.add_argument("--tmp-dir", default=None, help="Temporary directory root for CheckM2 working files (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
        parser.add_argument("--keep-temp", action="store_true", help="Keep successful CheckM2 temporary directories")
        parser.add_argument("--skip-checkm2", action="store_true", help="Skip CheckM2 during final bin refinement")
        parser.add_argument("--skip-refinement", action="store_true", help="Skip Binning_refiner combinations")
        parser.add_argument("--skip-consolidation", action="store_true", help="Skip final consolidation across bin sets")
        parser.add_argument("--keep-ambiguous", action="store_true", help="Keep ambiguous contigs in all bins")
        parser.add_argument("--remove-ambiguous", action="store_true", help="Remove ambiguous contigs from all bins")
        args = parser.parse_args()

        metacc_folder = args.metacc_folder
        bin3c_folder = args.bin3c_folder
        imputecc_folder = args.imputecc_folder
        metahict_folder = args.metahict_folder
        
        logger.info('Starting bin refinement with metahict')
        
        metahict = os.path.join(script_directory, 'bin_refinement.py')
        
        min_completeness = f"{args.min_completeness:g}"
        max_contamination = f"{args.max_contamination:g}"

        metahictCmd = [
            sys.executable,
            metahict,
            "-t", str(args.threads),
            "-o", metahict_folder,
            "-A", os.path.join(metacc_folder, "BIN"),
            "-B", os.path.join(imputecc_folder, "FINAL_BIN"),
            "-C", os.path.join(bin3c_folder, "fasta"),
            "-c", min_completeness,
            "-x", max_contamination,
            "--contamination-penalty", str(args.contamination_penalty),
            "--min-input-bin-size", str(args.min_input_bin_size),
            "--max-input-bin-size", str(args.max_input_bin_size),
            "--binning-refiner-min-size", str(args.binning_refiner_min_size),
        ]
        if args.tmp_dir:
            metahictCmd.extend(["--tmp-dir", args.tmp_dir])
        if args.keep_temp:
            metahictCmd.append("--keep-temp")
        if args.skip_checkm2:
            metahictCmd.append("--skip-checkm2")
        if args.skip_refinement:
            metahictCmd.append("--skip-refinement")
        if args.skip_consolidation:
            metahictCmd.append("--skip-consolidation")
        if args.keep_ambiguous:
            metahictCmd.append("--keep-ambiguous")
        if args.remove_ambiguous:
            metahictCmd.append("--remove-ambiguous")

        if os.path.isdir(metahict_folder):
            logger.info("Removing previous METAHICT integration output: %s", metahict_folder)
            shutil.rmtree(metahict_folder)
        os.makedirs(metahict_folder, exist_ok=True)

        final_bin_dir = os.path.join(metahict_folder, "final_bins")
        logger.info("metahictCmd : " + " ".join(metahictCmd))
        try:
            result = subprocess.run(
                metahictCmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True,
            )
            logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.error(e.stdout)
            raise
        
        logger.info('Bin refinement completed successfully')
        
        # Paths
        final_fasta = os.path.join(metahict_folder, "combined_final_bins.fa")
        contig_to_bin_tsv = os.path.join(metahict_folder, "contig_to_bin.tsv")
        
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
        
        # Write a portable, human-readable contig-to-bin table.
        with open(contig_to_bin_tsv, 'w') as f:
            f.write("contig\tbin_id\n")
            for contig, bin_id in final_clustering.items():
                f.write(f"{contig}\t{bin_id}\n")
        
        # The heatmap implementation requires a serialized mapping. Keep it
        # outside the public result set unless intermediate retention is
        # explicitly enabled by the caller.
        intermediate_folder = os.path.join(metahict_folder, "intermediates")
        os.makedirs(intermediate_folder, exist_ok=True)
        save_object(
            os.path.join(intermediate_folder, "final_bins.p"),
            final_clustering,
        )
        
    except Exception as e:
        logger.error(f'Error during bin refinement: {str(e)}')
        sys.exit(1)
