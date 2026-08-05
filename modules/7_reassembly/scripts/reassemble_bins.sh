#!/usr/bin/env bash

METAHICT_REASSEMBLE_ARGS=("$@")
module load anaconda3 2>/dev/null || true
CONDA_BASE="${CONDA_BASE:-/opt/anaconda3/2025.06-1}"
if [[ -s "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
	source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi
set -- "${METAHICT_REASSEMBLE_ARGS[@]}"

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
	echo "Usage: metahict reassemble_bins [options] -o output_dir -b bin_folder -1 reads_1.fastq -2 reads_2.fastq"
	echo ""
	echo "Options:"
	echo ""
	echo "	-b STR		folder with metagenomic bins"
	echo "	-o STR		output directory"
	echo "	-1 STR          forward reads to use for reassembly"
	echo "	-2 STR          reverse reads to use for reassembly"
	echo ""
	echo "	-t INT		number of threads (default=80)"
	echo "	-m INT		memory in GB (default=80% of available system memory)"
	echo "	-c INT		minimum desired bin completion % for final selection (default=50)"
	echo "	-x INT		maximum desired bin contamination % for final selection (default=10)"
	echo "	-l INT		minimum contig length to be included in reassembly (default=500)"
	echo ""
	echo "	--strict-cut-off	maximum allowed SNPs for strict read mapping (default=2)"
	echo "	--permissive-cut-off	maximum allowed SNPs for permissive read mapping (default=5)"
	echo "	--contamination-penalty FLOAT penalty used in final bin scoring (default=5)"
	echo "	--skip-checkm		do not run CheckM2 to assess bins"
	echo "	--parallel		run spades reassembly in parallel, but only using 1 thread per bin"
	echo "	--tmp-dir STR		temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)"
	echo "	--spades-mode STR	SPAdes mode: careful or none (default=careful)"
	echo "	--spades-phred-offset INT SPAdes phred offset; empty preserves SPAdes default"
	echo "	--spades-extra-args STR extra options passed to SPAdes (default=empty)"
	echo "	--skip-residual-assembly do not assemble unmapped reads with MEGAHIT"
	echo "	--keep-temp		keep successful SPAdes and CheckM2 temporary directories"
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

threads=80; mem=""; comp=50; cont=10; len=500
bins=None; f_reads=None; r_reads=None; out=None
strict_max=2; permissive_max=5
contamination_penalty=5
run_checkm=true
run_parallel=false
checkm2_db=""
tmp_root="${METAHICT_TMP_ROOT:-${TMPDIR:-/tmp}}"
if [[ ! -d "$tmp_root" || ! -w "$tmp_root" ]]; then
	tmp_root="/tmp"
fi
if [[ ! -d "$tmp_root" || ! -w "$tmp_root" ]]; then
	error "No writable temporary directory found. Set METAHICT_TMP_ROOT or TMPDIR to a writable short path."
fi
spades_phred_offset="${METAHICT_SPADES_PHRED_OFFSET:-}"
spades_mode="careful"
spades_extra_args=""
run_residual_assembly=true
keep_temp=false

while [[ $# -gt 0 ]]; do
	case "$1" in
		-t) threads=$2; shift 2;;
		-m) mem=$2; shift 2;;
		-o) out=$2; shift 2;;
		-b) bins=$2; shift 2;;
		-c) comp=$2; shift 2;;
		-x) cont=$2; shift 2;;
		-l) len=$2; shift 2;;
		-1) f_reads=$2; shift 2;;
		-2) r_reads=$2; shift 2;;
		-h|--help) help_message; exit 0;;
		--strict-cut-off) strict_max=$2; shift 2;;
		--permissive-cut-off) permissive_max=$2; shift 2;;
		--contamination-penalty) contamination_penalty=$2; shift 2;;
		--skip-checkm) run_checkm=false; shift 1;;
		--parallel) run_parallel=true; shift 1;;
		--checkm2_db) checkm2_db=$2; shift 2;;
		--tmp-dir) tmp_root=$2; shift 2;;
		--spades-mode) spades_mode=$2; shift 2;;
		--spades-phred-offset) spades_phred_offset=$2; shift 2;;
		--spades-extra-args) spades_extra_args=$2; shift 2;;
		--spades-extra-args=*) spades_extra_args="${1#*=}"; shift 1;;
		--skip-residual-assembly) run_residual_assembly=false; shift 1;;
		--keep-temp) keep_temp=true; shift 1;;
		*) echo "Unknown argument: $1"; help_message; exit 1;;
	esac
done

if [[ -z "$mem" ]]; then
	mem=$(awk '/MemTotal/ {printf "%d", ($2 / 1024 / 1024) * 0.8}' /proc/meminfo)
fi
if [[ "$spades_mode" != "careful" && "$spades_mode" != "none" ]]; then
	error "--spades-mode must be 'careful' or 'none'."
fi
if [[ ! -d "$tmp_root" || ! -w "$tmp_root" ]]; then
	tmp_root="/tmp"
fi
if [[ ! -d "$tmp_root" || ! -w "$tmp_root" ]]; then
	error "No writable temporary directory found. Set METAHICT_TMP_ROOT, TMPDIR, or --tmp-dir to a writable short path."
fi

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

if [ $out = None ] || [  $bins = None ] || [ $f_reads = None ] || [ $r_reads = None ] ; then 
	comm "Some non-optional parameters were not entered"
	help_message; exit 1
fi

for helper in print_comment.py filter_reads_for_bin_reassembly.py rm_short_contigs.py summarize_checkm2.py choose_best_bin.py plot_checkm2_results.py plot_reassembly.py; do
	if [ ! -s "${SOFT}/${helper}" ]; then
		error "Required helper script is missing: ${SOFT}/${helper}"
	fi
done

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
cp -r $bins ${out}/original_bins

if [ ! -d ${out}/binned_assembly ]; then mkdir ${out}/binned_assembly; fi

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
spades_status_dir=$(mktemp -d "${tmp_root%/}/metahict_spades_status.XXXXXX")

assemble () {
	bin_name=${1%_*}
	if [[ -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]]; then
		comm "Looks like $bin_name was already re-assembled. Skipping..."
	else
		tmp_dir=$(mktemp -d "${tmp_root%/}/metahict_spades_${bin_name}.XXXXXX")
		if [[ -z "$tmp_dir" || ! -d "$tmp_dir" ]]; then
			error "Could not create a short SPAdes temporary directory under ${tmp_root}."
		fi
		comm "NOW REASSEMBLING ${bin_name}"
		spades_cmd=(spades.py -t "$2" -m "$mem" --tmp "$tmp_dir")
		if [[ "$spades_mode" == "careful" ]]; then
			spades_cmd+=(--careful)
		fi
		if [[ -n "$spades_phred_offset" ]]; then
			spades_cmd+=(--phred-offset "$spades_phred_offset")
		fi
		if [[ -n "$spades_extra_args" ]]; then
			read -r -a spades_extra_args_array <<< "$spades_extra_args"
			spades_cmd+=("${spades_extra_args_array[@]}")
		fi
		spades_cmd+=(--untrusted-contigs "${out}/original_bins/${bin_name%.*}.fa")
		spades_cmd+=(-1 "${out}/reads_for_reassembly/${1%_*}_1.fastq")
		spades_cmd+=(-2 "${out}/reads_for_reassembly/${1%_*}_2.fastq")
		spades_cmd+=(-o "${out}/reassemblies/${bin_name}")
		"${spades_cmd[@]}"
		spades_status=$?
		if [[ $spades_status -eq 0 && "$keep_temp" = false ]]; then
			rm -rf "$tmp_dir"
		else
			if [[ $spades_status -eq 0 ]]; then
				comm "SPAdes temporary files kept at ${tmp_dir}"
			else
				warning "SPAdes failed for ${bin_name}; temporary files kept at ${tmp_dir}"
				touch "${spades_status_dir}/${bin_name}.failed"
				return $spades_status
			fi
		fi
	fi
}

if [ "$run_parallel" = true ]; then
	open_sem $threads
	for i in $(ls ${out}/reads_for_reassembly/ | grep _1.fastq); do 
		# The semaphore limits the number of simultaneous bin reassemblies to
		# the requested thread count, so each parallel SPAdes job must receive
		# one thread.  Passing $threads here would oversubscribe the allocation
		# by launching up to threads jobs that each request threads CPUs.
		run_with_lock assemble "$i" 1
	done
	wait
	if [[ $(find "$spades_status_dir" -type f -name '*.failed' | wc -l) -gt 0 ]]; then
		find "$spades_status_dir" -type f -name '*.failed' -printf '%f\n' >&2
		rm -rf "$spades_status_dir"
		error "One or more SPAdes reassemblies failed."
	fi
	comm "all assemblies complete"
else
	for i in $(ls ${out}/reads_for_reassembly/ | grep _1.fastq); do
		if ! assemble $i $threads; then
			rm -rf "$spades_status_dir"
			error "SPAdes reassembly failed."
		fi
	done
	comm "all assemblies complete"
fi
rm -rf "$spades_status_dir"

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
########################       ASSEMBLE RESIDUAL READS NOT RECRUITED TO BINS       ####################
########################################################################################################
announcement "ASSEMBLE RESIDUAL READS"

residual_contigs="${out}/residual_contigs.fa"
residual_assembly_dir="${out}/residual_assembly"
rm -f "$residual_contigs"

if [ "$run_residual_assembly" = false ]; then
	warning "Residual-read assembly skipped by --skip-residual-assembly."
	touch "$residual_contigs"
elif [[ ! -s "${out}/unmapped_shotgun_1.fastq" || ! -s "${out}/unmapped_shotgun_2.fastq" ]]; then
	warning "No unmapped read pairs available for residual assembly."
	touch "$residual_contigs"
else
	if ! command -v megahit >/dev/null 2>&1; then
		error "MEGAHIT is required for residual-read assembly but was not found in PATH."
	fi
	rm -rf "$residual_assembly_dir"
	megahit_tmp=$(mktemp -d "${tmp_root%/}/metahict_megahit_residual.XXXXXX")
	if [[ -z "$megahit_tmp" || ! -d "$megahit_tmp" ]]; then
		error "Could not create a short MEGAHIT temporary directory under ${tmp_root}."
	fi
	if megahit \
		-1 "${out}/unmapped_shotgun_1.fastq" \
		-2 "${out}/unmapped_shotgun_2.fastq" \
		-o "$residual_assembly_dir" \
		--tmp-dir "$megahit_tmp" \
		--num-cpu-threads "$threads" \
		--min-contig-len "$len" \
		--k-min 21 \
		--k-max 141 \
		--k-step 12 \
		--merge-level 20,0.95
	then
		megahit_status=0
	else
		megahit_status=$?
	fi
	if [[ "$keep_temp" = false ]]; then
		rm -rf "$megahit_tmp"
	else
		comm "MEGAHIT temporary files kept at $megahit_tmp"
	fi
	if [[ $megahit_status -ne 0 ]]; then
		error "MEGAHIT residual-read assembly failed."
	fi
	if [[ -s "${residual_assembly_dir}/final.contigs.fa" ]]; then
		cp "${residual_assembly_dir}/final.contigs.fa" "$residual_contigs"
		comm "Residual contigs saved to ${residual_contigs}"
	else
		warning "MEGAHIT completed but produced no residual contigs."
		touch "$residual_contigs"
	fi
fi

########################################################################################################
########################           RUN CHECKM2 (OPTIONAL CUSTOM DB)            ########################
########################################################################################################
if [ "$run_checkm" = true ]; then
	announcement "RUN CHECKM2 ON REASSEMBLED AND ORIGINAL BIN CANDIDATES"

	if [[ -n "$checkm2_db" ]]; then
		export CHECKM2DB="$checkm2_db"
		comm "Using custom CheckM2 database: $CHECKM2DB"
	else
		comm "Using default CheckM2 database"
	fi

	for base in $(ls ${out}/original_bins/ | grep "\.fa$"); do
		cp "${out}/original_bins/${base}" "${out}/reassembled_bins/${base%.*}.orig.fa"
	done

	if [[ -d ${out}/reassembled_bins.checkm2 ]]; then rm -r ${out}/reassembled_bins.checkm2; fi
	checkm2_tmp=$(mktemp -d "${tmp_root%/}/metahict_checkm2.XXXXXX")
	if [[ -z "$checkm2_tmp" || ! -d "$checkm2_tmp" ]]; then
		error "Could not create a short CheckM2 temporary directory under ${tmp_root}."
	fi
	old_tmpdir="${TMPDIR:-}"
	export TMPDIR="$checkm2_tmp"
	conda deactivate
	conda activate checkm2
	checkm2 predict -i ${out}/reassembled_bins -o ${out}/reassembled_bins.checkm2 -x fa -t $threads --tmpdir "$checkm2_tmp"
	checkm2_status=$?
	conda deactivate
	conda activate metahict_env
	if [[ -n "$old_tmpdir" ]]; then
		export TMPDIR="$old_tmpdir"
	else
		unset TMPDIR
	fi
	if [[ "$keep_temp" = false ]]; then
		rm -rf "$checkm2_tmp"
	else
		comm "CheckM2 temporary files kept at $checkm2_tmp"
	fi
	if [[ $checkm2_status -ne 0 || ! -s ${out}/reassembled_bins.checkm2/quality_report.tsv ]]; then error "CheckM2 failed."; fi

	${SOFT}/summarize_checkm2.py "${out}/reassembled_bins.checkm2/quality_report.tsv" | (read -r; printf "%s\n" "$REPLY"; sort) > "${out}/reassembled_bins.stats"
	if [[ ! -s ${out}/reassembled_bins.stats ]]; then error "Cannot make CheckM2 summary file."; fi

	########################################################################################################
	########################          FINDING THE BEST VERSION OF EACH BIN          ########################
	########################################################################################################
	announcement "FINDING THE BEST VERSION OF EACH BIN"

	rm -rf "${out}/reassembled_best_bins"
	mkdir "${out}/reassembled_best_bins"
	name_map="${out}/reassembled_bin_name_map.tsv"
	printf "final_bin\tselected_candidate\tselection_type\n" > "$name_map"
	selected_orig_count=0
	selected_strict_count=0
	selected_permissive_count=0
	for i in $(${SOFT}/choose_best_bin.py "${out}/reassembled_bins.stats" "$comp" "$cont" "$contamination_penalty"); do
		if [[ -s "${out}/reassembled_bins/${i}.fa" ]]; then
			selected_type="${i##*.}"
			base_bin="${i%.*}"
			final_bin="${base_bin//./}"
			case "$selected_type" in
				orig) selected_orig_count=$((selected_orig_count + 1)) ;;
				strict) selected_strict_count=$((selected_strict_count + 1)) ;;
				permissive) selected_permissive_count=$((selected_permissive_count + 1)) ;;
			esac
			cp "${out}/reassembled_bins/${i}.fa" "${out}/reassembled_best_bins/${final_bin}.fa"
			printf "%s\t%s\t%s\n" "$final_bin" "$i" "$selected_type" >> "$name_map"
		else
			warning "Best-bin candidate was selected but missing: ${out}/reassembled_bins/${i}.fa"
		fi
	done

	o=$selected_orig_count
	s=$selected_strict_count
	p=$selected_permissive_count
	announcement "Reassembly results are in! $s bins were improved with strict reassembly, $p bins were improved with permissive reassembly, and $o bins were not improved by any reassembly."

	if [[ $(ls "${out}/reassembled_best_bins" | wc -l) -lt 1 ]]; then
		error "There are no good bins in ${out}/reassembled_best_bins."
	fi

	rm -rf "${out}/work_files"
	mkdir "${out}/work_files"
	mv "${out}/reassembled_bins" "${out}/work_files/"
	mv "${out}/reassembled_bins.stats" "${out}/work_files/"
	cp "$name_map" "${out}/work_files/reassembled_bin_name_map.tsv"
	for work_dir in reads_for_reassembly binned_assembly reassemblies; do
		if [[ -e "${out}/${work_dir}" ]]; then
			mv "${out}/${work_dir}" "${out}/work_files/"
		fi
	done
	mv "${out}/reassembled_best_bins" "${out}/reassembled_bins"

	announcement "RE-RUN CHECKM2 ON FINAL BEST REASSEMBLED BINS"
	if [[ -d ${out}/reassembled_bins.checkm2 ]]; then rm -r ${out}/reassembled_bins.checkm2; fi
	checkm2_tmp=$(mktemp -d "${tmp_root%/}/metahict_checkm2.XXXXXX")
	if [[ -z "$checkm2_tmp" || ! -d "$checkm2_tmp" ]]; then
		error "Could not create a short CheckM2 temporary directory under ${tmp_root}."
	fi
	old_tmpdir="${TMPDIR:-}"
	export TMPDIR="$checkm2_tmp"
	conda deactivate
	conda activate checkm2
	checkm2 predict -i "${out}/reassembled_bins" -o "${out}/reassembled_bins.checkm2" -x fa -t "$threads" --tmpdir "$checkm2_tmp"
	checkm2_status=$?
	conda deactivate
	conda activate metahict_env
	if [[ -n "$old_tmpdir" ]]; then
		export TMPDIR="$old_tmpdir"
	else
		unset TMPDIR
	fi
	if [[ "$keep_temp" = false ]]; then
		rm -rf "$checkm2_tmp"
	else
		comm "CheckM2 temporary files kept at $checkm2_tmp"
	fi
	if [[ $checkm2_status -ne 0 || ! -s ${out}/reassembled_bins.checkm2/quality_report.tsv ]]; then error "CheckM2 failed."; fi

	${SOFT}/summarize_checkm2.py "${out}/reassembled_bins.checkm2/quality_report.tsv" | (read -r; printf "%s\n" "$REPLY"; sort) > "${out}/reassembled_bins.stats"
	if [[ ! -s ${out}/reassembled_bins.stats ]]; then error "Cannot make final CheckM2 summary file."; fi
	${SOFT}/plot_checkm2_results.py "${out}/reassembled_bins.checkm2/quality_report.tsv" "${out}/reassembled_bins"
	head -n 1 "${out}/work_files/reassembled_bins.stats" > "${out}/original_bins.stats"
	grep "orig" "${out}/work_files/reassembled_bins.stats" >> "${out}/original_bins.stats" || true
	${SOFT}/plot_reassembly.py "$out" "$comp" "$cont" "${out}/reassembled_bins.stats" "${out}/original_bins.stats"
else
	warning "CheckM2 was skipped, so final best-bin selection was not performed."
fi

comm "Writing combined contigs FASTA"
export REASSEMBLED_BINS_DIR="${out}/reassembled_bins"
export RESIDUAL_CONTIGS_FASTA="${residual_contigs}"
export COMBINED_CONTIGS_FASTA="${out}/combined_contigs.fa"
python3 <<'EOF'
import os

bins_dir = os.environ["REASSEMBLED_BINS_DIR"]
residual_fasta = os.environ["RESIDUAL_CONTIGS_FASTA"]
combined_fasta = os.environ["COMBINED_CONTIGS_FASTA"]
seen = {}

def read_fasta(path):
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, chunks
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, chunks

with open(combined_fasta, "w") as out_f:
    for filename in sorted(os.listdir(bins_dir)):
        if not filename.endswith((".fa", ".fasta", ".fna")):
            continue
        bin_name = os.path.splitext(filename)[0]
        fasta = os.path.join(bins_dir, filename)
        for old_id, seq_chunks in read_fasta(fasta):
            new_id = f"{bin_name}|{old_id}"
            copy_number = seen.get(new_id, 0) + 1
            seen[new_id] = copy_number
            if copy_number > 1:
                new_id = f"{new_id}|copy{copy_number}"
            out_f.write(f">{new_id}\n")
            for chunk in seq_chunks:
                out_f.write(f"{chunk}\n")
    if residual_fasta and os.path.exists(residual_fasta):
        for old_id, seq_chunks in read_fasta(residual_fasta):
            new_id = f"residual|{old_id}"
            copy_number = seen.get(new_id, 0) + 1
            seen[new_id] = copy_number
            if copy_number > 1:
                new_id = f"{new_id}|copy{copy_number}"
            out_f.write(f">{new_id}\n")
            for chunk in seq_chunks:
                out_f.write(f"{chunk}\n")
EOF
if [[ ! -s ${out}/combined_contigs.fa ]]; then
	error "Combined contigs FASTA was not created."
fi

########################################################################################################
########################          REPORT CIRCULAR REASSEMBLED CONTIGS          ########################
########################################################################################################
announcement "REPORT CIRCULAR REASSEMBLED CONTIGS"

export REASSEMBLED_BINS_DIR="${out}/reassembled_bins"
export CIRCULAR_CONTIGS_REPORT="${out}/circular_contigs.tsv"
python3 <<'EOF'
import os

bins_dir = os.environ["REASSEMBLED_BINS_DIR"]
report = os.environ["CIRCULAR_CONTIGS_REPORT"]
min_overlap = 50
max_overlap = 1000

def read_fasta(path):
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks).upper()

def terminal_overlap(seq):
    max_len = min(max_overlap, len(seq) // 2)
    for size in range(max_len, min_overlap - 1, -1):
        if seq[:size] == seq[-size:]:
            return size
    return 0

rows = []
for fname in sorted(os.listdir(bins_dir)):
    if not fname.endswith((".fa", ".fasta", ".fna")):
        continue
    fasta = os.path.join(bins_dir, fname)
    bin_id = os.path.splitext(fname)[0]
    for contig_id, seq in read_fasta(fasta):
        overlap = terminal_overlap(seq)
        rows.append((
            contig_id,
            "reassembled_bin",
            bin_id,
            str(len(seq)),
            "circular" if overlap else "not_circular",
            "terminal_overlap" if overlap else "none",
            str(overlap),
        ))

with open(report, "w") as out:
    out.write("contig_id\tsource\tbin_id\tlength\tcircularity_status\tevidence\toverlap_length\n")
    for row in rows:
        out.write("\t".join(row) + "\n")

total = len(rows)
circular = sum(1 for row in rows if row[4] == "circular")
print(f"[INFO] Circular contig report saved to {report}")
print(f"[INFO] Circular contigs detected: {circular} of {total}")
EOF

########################################################################################################
########################    REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!        ########################
########################################################################################################
announcement "BIN REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!"
