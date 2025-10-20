#!/usr/bin/env bash

module load anaconda3

# Set the path to helper scripts relative to the script location
SOFT="$(dirname "$(realpath "$0")")/bin_integration"

##############################################################################################################################################################
#
# Improves a set of bins by aligning reads back to the bins and reassembling them.
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################

help_message () {
	echo ""
	echo "Usage: metahit reassemble_bins [options] -o output_dir -b bin_folder -1 reads_1.fastq -2 reads_2.fastq"
	echo ""
	echo "Options:"
	echo ""
	echo "	-b STR		folder with metagenomic bins"
	echo "	-o STR		output directory"
	echo "	-1 STR          forward reads to use for reassembly"
	echo "	-2 STR          reverse reads to use for reassembly"
	echo ""
	echo "	--nanopore STR	nanopore reads to use for reassembly"
	echo ""
	echo "	-t INT		number of threads (default=1)"
	echo "	-m INT		memory in GB (default=120)"
	echo "	-c INT		minimum desired bin completion % (default=70)"
	echo "	-x INT		maximum desired bin contamination % (default=10)"
	echo "	-l INT		minimum contig length to be included in reassembly (default=500)"
	echo ""
	echo "	--strict-cut-off	maximum allowed SNPs for strict read mapping (default=2)"
	echo "	--permissive-cut-off	maximum allowed SNPs for permissive read mapping (default=5)"
	echo "	--skip-checkm		dont run CheckM to assess bins"
	echo "	--parallel		run spades reassembly in parallel, but only using 1 thread per bin"
	echo "	--mdmcleaner		the bin directory have results from MDMcleaner"
	echo "	--checkm2_db STR	path to CheckM2 database (optional)"
	echo ""
}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for i in $(seq $1 -1 1); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    "$@" 
    printf '%.3d' $? >&3
    )&
}

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

threads=1; mem=120; comp=50; cont=10; len=500
bins=None; f_reads=None; r_reads=None; out=None
strict_max=2; permissive_max=5
run_checkm=true
run_parallel=false
nanopore=false
mdmcleaner=false
checkm2_db=""

OPTS=`getopt -o ht:m:o:x:c:l:b:1:2: --long help,parallel,skip-checkm,strict-cut-off:,permissive-cut-off:,nanopore:,mdmcleaner,checkm2_db: -- "$@"`
if [ $? -ne 0 ]; then help_message; exit 1; fi

while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
                -m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
		-x) cont=$2; shift 2;;
		-c) comp=$2; shift 2;;
		-b) bins=$2; shift 2;;
		-l) len=$2; shift 2;;
		-1) f_reads=$2; shift 2;;
		-2) r_reads=$2; shift 2;;
                -h | --help) help_message; exit 0; shift 1;;
		--strict-cut-off) strict_max=$2; shift 2;;
		--permissive-cut-off) permissive_max=$2; shift 2;;
		--skip-checkm) run_checkm=false; shift 1;;
		--parallel) run_parallel=true; shift 1;;
		--nanopore) nanopore_reads=$2;nanopore=true; shift 2;;
		--mdmcleaner) mdmcleaner=true; shift 1;;
		--checkm2_db) checkm2_db=$2; shift 2;;
                --) shift; break ;;
                *) break;;
        esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

if [ $out = None ] || [  $bins = None ] || [ $f_reads = None ] || [ $r_reads = None ] ; then 
	comm "Some non-optional parameters were not entered"
	help_message; exit 1
fi

if [ ! -s ${SOFT}/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same folder as the main scripts."
fi

########################################################################################################
########################               BEGIN REASSEMBLY PIPELINE!               ########################
########################################################################################################
announcement "BEGIN PIPELINE!"
comm "setting up output folder and copying over bins..."
if [ ! -d $out ]; then
        mkdir $out;
else
        warning "Warning: $out already exists!"
fi

if [ -d ${out}/original_bins ]; then rm -r ${out}/original_bins; fi
if [ "$mdmcleaner" = true ]; then
	mkdir ${out}/original_bins
	for i in $bins/*/*_kept_contigs.fasta.gz; do gunzip -k $i && mv ${i%.gz} ${out}/original_bins; done
else
	cp -r $bins ${out}/original_bins
fi

if [ ! -d ${out}/binned_assembly ]; then mkdir ${out}/binned_assembly; fi

if [ "$mdmcleaner" = true ]; then
comm "Cleaning mdmcleaner contigs..."
for i in ${out}/original_bins/*_filtered_kept_contigs.fasta
	do  awk -F " " '/^>/ { print $1; next } 1' $i > tmp.file && \
	mv tmp.file ${i%_filtered_kept_contigs.fasta}.fa && \
	rm $i
done
rm -f tmp.file
fi

if [ -s ${out}/binned_assembly/assembly.fa ]; then rm ${out}/binned_assembly/assembly.fa; fi
for i in $(ls ${out}/original_bins); do cat ${out}/original_bins/$i >> ${out}/binned_assembly/assembly.fa; done

########################################################################################################
########################        RECRUITING READS TO BINS FOR REASSEMBLY         ########################
########################################################################################################
announcement "RECRUITING READS TO BINS FOR REASSEMBLY"

ulimit -n 10000
if [[ $? -ne 0 ]]; then
	ULIMIT=$(ulimit -n)
	warning "Your OS allows $ULIMIT open files. You may need fewer bins."
fi

if [[ ! -s ${out}/binned_assembly/assembly.fa.amb ]]; then
	comm "Indexing the assembly"
	bwa index ${out}/binned_assembly/assembly.fa
	if [[ $? -ne 0 ]]; then error "BWA failed to index $i"; fi

	if [ -d ${out}/reads_for_reassembly ]; then rm -r ${out}/reads_for_reassembly; fi
	mkdir ${out}/reads_for_reassembly

	comm "Aligning all reads back to entire assembly, extracting unmapped reads, and splitting by bin"

    bwa mem -t $threads ${out}/binned_assembly/assembly.fa $f_reads $r_reads \
    | tee >(samtools view -b -f 12 | samtools fastq -1 ${out}/unmapped_shotgun_1.fastq -2 ${out}/unmapped_shotgun_2.fastq -0 /dev/null -s /dev/null) \
    | ${SOFT}/filter_reads_for_bin_reassembly.py ${out}/original_bins ${out}/reads_for_reassembly $strict_max $permissive_max
    
    comm "Unmapped reads saved to ${out}/unmapped_shotgun_1.fastq and ${out}/unmapped_shotgun_2.fastq"
else
	comm "Assembly was already indexed. Skipping indexing and splitting."
fi

########################################################################################################
########################             REASSEMBLING BINS WITH SPADES              ########################
########################################################################################################
announcement "REASSEMBLING BINS WITH SPADES"
mkdir ${out}/reassemblies

assemble () {
	bin_name=${1%_*}
	if [[ -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]]; then
		comm "Looks like $bin_name was already re-assembled. Skipping..."
	else
		tmp_dir=${out}/reassemblies/${bin_name}.tmp
		mkdir $tmp_dir
		comm "NOW REASSEMBLING ${bin_name}"
		spades.py -t $2 -m $mem --tmp $tmp_dir --careful \
		--untrusted-contigs ${out}/original_bins/${bin_name%.*}.fa \
		-1 ${out}/reads_for_reassembly/${1%_*}_1.fastq \
		-2 ${out}/reads_for_reassembly/${1%_*}_2.fastq \
		-o ${out}/reassemblies/${bin_name}
	fi
}

if [ "$run_parallel" = true ]; then
	open_sem $threads
	for i in $(ls ${out}/reads_for_reassembly/ | grep _1.fastq); do 
		run_with_lock assemble $i $threads 
	done
	wait
	comm "all assemblies complete"
else
	for i in $(ls ${out}/reads_for_reassembly/ | grep _1.fastq); do
		assemble $i $threads 
	done
	comm "all assemblies complete"
fi

comm "Finalizing reassemblies"
mkdir ${out}/reassembled_bins
for i in $( ls ${out}/reassemblies/ ); do
	spades_folder=${out}/reassemblies/$i
	bin_name=${spades_folder##*/}
	if [ -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]; then
		${SOFT}/rm_short_contigs.py $len ${out}/reassemblies/${bin_name}/scaffolds.fasta > ${out}/reassembled_bins/${bin_name}.fa
	fi
done

if [[ $(ls ${out}/reassembled_bins/ | wc -l) -lt 1 ]]; then
	error "None of the bins were successfully reassembled."
fi

########################################################################################################
########################           RUN CHECKM2 (OPTIONAL CUSTOM DB)            ########################
########################################################################################################
if [ "$run_checkm" = true ]; then
	announcement "RUN CHECKM2 ON REASSEMBLED BINS"

	if [[ -n "$checkm2_db" ]]; then
		export CHECKM2DB="$checkm2_db"
		comm "Using custom CheckM2 database: $CHECKM2DB"
	else
		comm "Using default CheckM2 database"
	fi

	if [[ -d ${out}/reassembled_bins.checkm2 ]]; then rm -r ${out}/reassembled_bins.checkm2; fi
	mkdir ${out}/tmp
	conda deactivate
	conda activate checkm2
	checkm2 predict -i ${out}/reassembled_bins -o ${out}/reassembled_bins.checkm2 -x fa -t $threads --tmpdir ${out}/tmp
	conda deactivate
	conda activate metahit_env
	if [[ ! -s ${out}/reassembled_bins.checkm2/quality_report.tsv ]]; then error "CheckM2 failed."; fi
	rm -r ${out}/tmp
fi

########################################################################################################
########################    REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!        ########################
########################################################################################################
announcement "BIN REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!"
