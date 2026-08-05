#!/usr/bin/env python
import argparse
import subprocess
import os
import sys

script_dir = sys.path[0]

def run_command(command):
    print(f"[INFO] Executing command:\n{command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {e}")
        sys.exit(1)

def ensure_dir_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def absolute_path(path):
    return os.path.abspath(path)

def preprocessing_defaults_text():
    return """Module 1 preprocessing defaults:
threads=80
prefix=derived from input R1 filename
dedup=false
minlen=50
trimq=10
qtrim=r
ftl=10
xmx=80% of available system memory
ftm=5
ktrim=r
k=23
mink=11
hdist=1
adapter_ref=<project_path>/external/bbmap/resources/adapters.fa
skip_pre_qc_report=false
skip_post_qc_report=false"""

def assembly_defaults_text():
    return """Module 2 assembly defaults:
threads=80
memory=80% of available system memory
min_len=1000
assembler=megahit
megahit_k_min=21
megahit_k_max=141
megahit_k_step=12
megahit_merge_level=20,0.95
metaspades_k_list=21,33,55
flye_method=--nano-raw
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
skip_quast=false
keep_temp=false"""

def alignment_defaults_text():
    return """Module 3 alignment defaults:
threads=80
bwa_options=-5SP
samtools_filter=-F 0x900
mapq=30
min_intra_dist=10000
min_match_len=30
sort_memory=1G
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
keep_sam=false
skip_metrics=false"""

def coverage_defaults_text():
    return """Module 4 coverage defaults:
threads=80
memory=80% of available system memory
percent_identity=97
min_mapq=0
weight_mapq=0.0
include_edge_bases=false
max_edge_bases=75
min_contig_length=0
min_contig_depth=0
bbmap_extra_args=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
keep_sam=false
keep_temp=false"""

def contact_defaults_text():
    return """Module 5 contact defaults:
normalization_method=required
metacc_min_signal=1
metacc_min_len=1000
metacc_min_mapq=30
metacc_min_match=30
spurious_contact_percent=5
coverage_file=
epsilon=1
max_iter=1000
tol=1e-6"""

def binning_defaults_text():
    return """Module 6 binning defaults:
threads=80
enzyme=Sau3AI,MluCI
metacc_min_len=1000
metacc_min_signal=2
metacc_min_mapq=30
metacc_min_match=30
metacc_min_binsize=150000
normcc_thres=0.05
bin3c_min_len=1000
bin3c_min_signal=5
bin3c_min_mapq=60
bin3c_min_match=10
bin3c_min_extent=50000
min_completeness=50
max_contamination=10
contamination_penalty=5
min_input_bin_size=50000
max_input_bin_size=20000000
binning_refiner_min_size=524288
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
no_fasta=false
no_report=false
no_spades=false
only_large=false
keep_temp=false
skip_checkm2=false
skip_refinement=false
skip_consolidation=false
keep_ambiguous=false
remove_ambiguous=false
seed="""

def reassembly_defaults_text():
    return """Module 7 reassembly defaults:
threads=80
memory=80% of available system memory
cutoff_quantile=0.95
top_k=100
min_mapq=30
min_match_len=30
exclude_duplicates=false
write_nonselected_hic=false
min_contig_len=500
strict_cut_off=2
permissive_cut_off=5
contamination_penalty=5
skip_checkm2=false
checkm2_db=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
spades_mode=careful
spades_phred_offset=
spades_extra_args=
skip_residual_assembly=false
keep_temp=false"""

def scaffolding_defaults_text():
    return """Module 8 scaffolding defaults:
threads=80
memory=80% of available system memory
resolution=10000
min_contig_len=5000
bwa_options=-5SP
samtools_filter=-F 0x900
sort_memory=1G
metacc_min_mapq=30
metacc_min_len=1000
metacc_min_match=30
metacc_min_signal=2
bin3c_min_mapq=60
bin3c_min_len=1000
bin3c_min_match=10
bin3c_min_signal=5
yahs_resolutions=
yahs_min_mapq=10
yahs_min_contig_len=0
yahs_rounds=1
yahs_no_contig_ec=false
yahs_no_scaffold_ec=false
yahs_no_mem_check=false
yahs_extra_args=
normcc_thres=0.05
heatmap_max_image=5000
skip_checkm2=false
checkm2_db=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
keep_temp=false"""

def annotation_defaults_text():
    return """Module 9 annotation defaults:
threads=80
pplacer_cpus=same as threads
gtdbtk_db=<project_path>/databases/gtdbtk_db/release220
extension=fa
prefix=gtdbtk
skip_ani_screen=true
mash_db=
min_perc_aa=10
min_af=0.5
full_tree=false
scratch_dir=
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp
force=false
keep_intermediates=false
debug=false
write_single_copy_genes=false
gtdbtk_extra_args="""

def mge_defaults_text():
    return """Module 10 MGE defaults:
threads=80
genomad_db=<project_path>/databases/genomad_db/genomad_db
checkv_db=<project_path>/databases/checkv_db/checkv-db-v1.5
genomad_splits=8
genomad_sensitivity=4.2
genomad_cleanup=true
genomad_restart=false
genomad_preset=default
genomad_min_score=0.7
genomad_max_fdr=0.1
genomad_extra_args=
checkv_remove_tmp=false
checkv_restart=false
checkv_trimming=default
checkv_extra_args=
association_filter=zscore
zscore_threshold=0.5
fixed_contact_threshold=0
top_percent=50
min_raw_contacts=2
viral_quality_levels=Complete,High-quality,Medium-quality,Low-quality
min_terminal_overlap=50
max_terminal_overlap=1000
min_contact_strength=0
tmp_dir=METAHICT_TMP_ROOT, TMPDIR, or /tmp"""

def preprocessing(args):
    if args.print_defaults:
        print(preprocessing_defaults_text())
        return

    print("[INFO] Running Preprocessing Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/1_preprocessing/1_preprocessing.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )
    
    if args.prefix:
        command += f' --prefix "{args.prefix}"'
    command += (
        f' --minlen {args.minlen}'
        f' --trimq {args.trimq}'
        f' --qtrim "{args.qtrim}"'
        f' --ftl {args.ftl}'
        f' --ftm {args.ftm}'
        f' --ktrim "{args.ktrim}"'
        f' --k {args.k}'
        f' --mink {args.mink}'
        f' --hdist {args.hdist}'
    )
    if args.xmx:
        command += f' --xmx "{args.xmx}"'
    if args.adapter_ref:
        command += f' --adapter-ref "{absolute_path(args.adapter_ref)}"'
    if args.dedup:
        command += " --dedup"
    if args.skip_pre_qc_report:
        command += " --skip-pre-qc-report"
    if args.skip_post_qc_report:
        command += " --skip-post-qc-report"

    run_command(command)
    
def assembly(args):
    if args.print_defaults:
        print(assembly_defaults_text())
        return

    print("[INFO] Running Assembly Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/2_assembly/2_assembly.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.memory:
        command += f' -m {args.memory}'
    if args.min_len:
        command += f' -l {args.min_len}'
    if args.megahit:
        command += " --megahit"
    if args.metaspades:
        command += " --metaspades"
    if args.metaflye:
        command += " --metaflye"
    command += (
        f' --k-min {args.k_min}'
        f' --k-max {args.k_max}'
        f' --k-step {args.k_step}'
        f' --merge-level "{args.merge_level}"'
        f' --k-list "{args.k_list}"'
        f' --flye-method "{args.flye_method}"'
    )
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.skip_quast:
        command += " --skip-quast"
    if args.keep_temp:
        command += " --keep-temp"

    run_command(command)
    
def alignment(args):
    if args.print_defaults:
        print(alignment_defaults_text())
        return

    print("[INFO] Running Alignment Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/3_alignment/3_alignment.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-r "{absolute_path(args.reference)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.samtools_filter:
        command += f' --samtools-filter "{args.samtools_filter}"'
    command += (
        f' --bwa-options "{args.bwa_options}"'
        f' --mapq {args.mapq}'
        f' --min-intra-dist {args.min_intra_dist}'
        f' --min-match-len {args.min_match_len}'
        f' --sort-memory "{args.sort_memory}"'
    )
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.keep_sam:
        command += " --keep-sam"
    if args.skip_metrics:
        command += " --skip-metrics"

    run_command(command)
    
def coverage(args):
    if args.print_defaults:
        print(coverage_defaults_text())
        return

    print("[INFO] Running Coverage Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/4_coverage/4_coverage.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'-1 "{absolute_path(args.r1)}" '
        f'-2 "{absolute_path(args.r2)}" '
        f'-r "{absolute_path(args.reference)}" '
        f'-o "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.memory:
        command += f' -m "{args.memory}"'
    command += (
        f' --percent-identity {args.percent_identity}'
        f' --min-mapq {args.min_mapq}'
        f' --weight-mapq {args.weight_mapq}'
        f' --max-edge-bases {args.max_edge_bases}'
        f' --min-contig-length {args.min_contig_length}'
        f' --min-contig-depth {args.min_contig_depth}'
    )
    if args.include_edge_bases:
        command += " --include-edge-bases"
    if args.bbmap_extra_args:
        command += f' --bbmap-extra-args "{args.bbmap_extra_args}"'
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.keep_sam:
        command += " --keep-sam"
    if args.keep_temp:
        command += " --keep-temp"

    run_command(command)
    
def contact(args):
    if args.print_defaults:
        print(contact_defaults_text())
        return

    print("[INFO] Running Contact Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/5_contact/5_contact.sh" '
        f'{args.method} '
        f'-p "{absolute_path(args.project_path)}" '
        f'--bam "{absolute_path(args.bam)}" '
        f'--fasta "{absolute_path(args.fasta)}" '
        f'--out "{output_dir}" '
        f'--enzyme "{args.enzyme}"'
    )

    command += (
        f' --metacc-min-signal {args.metacc_min_signal}'
        f' --metacc-min-len {args.metacc_min_len}'
        f' --metacc-min-mapq {args.metacc_min_mapq}'
        f' --metacc-min-match {args.metacc_min_match}'
        f' --spurious-contact-percent {args.spurious_contact_percent}'
        f' --epsilon {args.epsilon}'
        f' --max-iter {args.max_iter}'
        f' --tol {args.tol}'
    )
    if args.coverage_file:
        command += f' --coverage-file "{absolute_path(args.coverage_file)}"'

    run_command(command)
    
def binning(args):
    if args.print_defaults:
        print(binning_defaults_text())
        return

    required = {
        "--fasta": args.fasta,
        "--bam": args.bam,
        "--output": args.output,
        "--project_path": args.project_path,
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        print(f"[ERROR] Missing required argument(s): {', '.join(missing)}")
        sys.exit(1)

    print("[INFO] Running Binning Module")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/6_binning/6_binning.sh" '
        f'"{absolute_path(args.fasta)}" '
        f'"{absolute_path(args.bam)}" '
        f'"{output_dir}" '
        f'"{absolute_path(args.project_path)}" '
        f'--threads {args.threads}'
    )

    if args.checkm_db:
        command += f' --checkm_db "{absolute_path(args.checkm_db)}"'
    command += (
        f' --enzyme "{args.enzyme}"'
        f' --metacc-min-len {args.metacc_min_len}'
        f' --metacc-min-signal {args.metacc_min_signal}'
        f' --metacc-min-mapq {args.metacc_min_mapq}'
        f' --metacc-min-match {args.metacc_min_match}'
        f' --metacc-min-binsize {args.metacc_min_binsize}'
        f' --normcc-thres {args.normcc_thres}'
        f' --bin3c-min-len {args.bin3c_min_len}'
        f' --bin3c-min-signal {args.bin3c_min_signal}'
        f' --bin3c-min-mapq {args.bin3c_min_mapq}'
        f' --bin3c-min-match {args.bin3c_min_match}'
        f' --bin3c-min-extent {args.bin3c_min_extent}'
        f' --min-completeness {args.min_completeness}'
        f' --max-contamination {args.max_contamination}'
        f' --contamination-penalty {args.contamination_penalty}'
        f' --min-input-bin-size {args.min_input_bin_size}'
        f' --max-input-bin-size {args.max_input_bin_size}'
        f' --binning-refiner-min-size {args.binning_refiner_min_size}'
    )
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.no_fasta:
        command += " --no-fasta"
    if args.no_report:
        command += " --no-report"
    if args.no_spades:
        command += " --no-spades"
    if args.only_large:
        command += " --only-large"
    if args.keep_temp:
        command += " --keep-temp"
    if args.skip_checkm2:
        command += " --skip-checkm2"
    if args.skip_refinement:
        command += " --skip-refinement"
    if args.skip_consolidation:
        command += " --skip-consolidation"
    if args.keep_ambiguous:
        command += " --keep-ambiguous"
    if args.remove_ambiguous:
        command += " --remove-ambiguous"
    if args.seed is not None:
        command += f" --seed {args.seed}"

    run_command(command)
    
def reassembly(args):
    if args.print_defaults:
        print(reassembly_defaults_text())
        return

    print("[INFO] Running Reassembly Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/7_reassembly/7_reassembly.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--bin "{absolute_path(args.bin)}" '
        f'--assembly "{absolute_path(args.assembly)}" '
        f'--hic1 "{absolute_path(args.hic1)}" '
        f'--hic2 "{absolute_path(args.hic2)}" '
        f'--sg1 "{absolute_path(args.sg1)}" '
        f'--sg2 "{absolute_path(args.sg2)}" '
        f'--bam "{absolute_path(args.bam)}" '
        f'--outdir "{output_dir}" '
        f'-t {args.threads} '
        f'--cutoff-quantile {args.cutoff_quantile} '
        f'--top-k {args.top_k} '
        f'--min-mapq {args.min_mapq} '
        f'--min-match-len {args.min_match_len} '
        f'--min-contig-len {args.min_contig_len} '
        f'--strict-cut-off {args.strict_cut_off} '
        f'--permissive-cut-off {args.permissive_cut_off} '
        f'--contamination-penalty {args.contamination_penalty} '
        f'--spades-mode "{args.spades_mode}"'
    )

    if args.memory:
        command += f' -m {args.memory}'
    if args.exclude_duplicates:
        command += " --exclude-duplicates"
    if args.write_nonselected_hic:
        command += " --write-nonselected-hic"
    if args.skip_checkm2:
        command += " --skip-checkm2"
    if args.checkm2_db:
        command += f' --checkm2_db "{absolute_path(args.checkm2_db)}"'
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.spades_phred_offset:
        command += f' --spades-phred-offset {args.spades_phred_offset}'
    if args.spades_extra_args:
        command += f' --spades-extra-args "{args.spades_extra_args}"'
    if args.skip_residual_assembly:
        command += " --skip-residual-assembly"
    if args.keep_temp:
        command += " --keep-temp"

    run_command(command)
    
def scaffolding(args):
    if args.print_defaults:
        print(scaffolding_defaults_text())
        return

    print("[INFO] Running Scaffolding Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/8_scaffolding/8_scaffolding.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--fasta "{absolute_path(args.fasta)}" '
        f'--enzyme "{args.enzyme}" '
        f'--outdir "{output_dir}" '
        f'--hic1 "{absolute_path(args.hic1)}" '
        f'--hic2 "{absolute_path(args.hic2)}" '
        f'-t {args.threads} '
        f'-r {args.resolution} '
        f'--min-contig-len {args.min_contig_len} '
        f'--bwa-options "{args.bwa_options}" '
        f'--samtools-filter "{args.samtools_filter}" '
        f'--sort-memory "{args.sort_memory}" '
        f'--metacc-min-mapq {args.metacc_min_mapq} '
        f'--metacc-min-len {args.metacc_min_len} '
        f'--metacc-min-match {args.metacc_min_match} '
        f'--metacc-min-signal {args.metacc_min_signal} '
        f'--bin3c-min-mapq {args.bin3c_min_mapq} '
        f'--bin3c-min-len {args.bin3c_min_len} '
        f'--bin3c-min-match {args.bin3c_min_match} '
        f'--bin3c-min-signal {args.bin3c_min_signal} '
        f'--yahs-min-mapq {args.yahs_min_mapq} '
        f'--yahs-min-contig-len {args.yahs_min_contig_len} '
        f'--yahs-rounds {args.yahs_rounds} '
        f'--normcc-thres {args.normcc_thres} '
        f'--heatmap-max-image {args.heatmap_max_image}'
    )

    if args.bam:
        command += f' --bam "{absolute_path(args.bam)}"'
    if args.memory:
        command += f' -m {args.memory}'    
    if args.yahs_resolutions:
        command += f' --yahs-resolutions "{args.yahs_resolutions}"'
    if args.yahs_no_contig_ec:
        command += " --yahs-no-contig-ec"
    if args.yahs_no_scaffold_ec:
        command += " --yahs-no-scaffold-ec"
    if args.yahs_no_mem_check:
        command += " --yahs-no-mem-check"
    if args.yahs_extra_args:
        command += f' --yahs-extra-args "{args.yahs_extra_args}"'
    if args.skip_checkm2:
        command += " --skip-checkm2"
    if args.checkm2_db:
        command += f' --checkm2_db "{absolute_path(args.checkm2_db)}"'
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.keep_temp:
        command += " --keep-temp"
        
    run_command(command)
    
def annotation(args):
    if args.print_defaults:
        print(annotation_defaults_text())
        return

    print("[INFO] Running Annotation Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/9_annotation/9_annotation.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--bin "{absolute_path(args.bin)}" '
        f'--outdir "{output_dir}" '
        f'-t {args.threads}'
    )
    
    if args.pplacer_cpus:
        command += f' --pplacer-cpus {args.pplacer_cpus}'
    if args.gtdbtk_db:
        command += f' --gtdbtk_db "{absolute_path(args.gtdbtk_db)}"'
    if args.extension:
        command += f' --extension "{args.extension}"'
    if args.prefix:
        command += f' --prefix "{args.prefix}"'
    if args.skip_ani_screen:
        command += " --skip-ani-screen"
    else:
        command += " --no-skip-ani-screen"
    if args.mash_db:
        command += f' --mash-db "{absolute_path(args.mash_db)}"'
    if args.min_perc_aa is not None:
        command += f' --min-perc-aa {args.min_perc_aa}'
    if args.min_af is not None:
        command += f' --min-af {args.min_af}'
    if args.full_tree:
        command += " --full-tree"
    if args.scratch_dir:
        command += f' --scratch-dir "{absolute_path(args.scratch_dir)}"'
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'
    if args.force:
        command += " --force"
    if args.keep_intermediates:
        command += " --keep-intermediates"
    if args.debug:
        command += " --debug"
    if args.write_single_copy_genes:
        command += " --write-single-copy-genes"
    if args.gtdbtk_extra_args:
        command += f' --gtdbtk-extra-args "{args.gtdbtk_extra_args}"'

    run_command(command)

def mge(args):
    if args.print_defaults:
        print(mge_defaults_text())
        return

    print("[INFO] Running MGE Module")
    output_dir = absolute_path(args.outdir)
    ensure_dir_exists(output_dir)

    command = (
        f'"{script_dir}/modules/10_MGE/10_MGE.sh" '
        f'-p "{absolute_path(args.project_path)}" '
        f'--combined "{absolute_path(args.combined)}" '
        f'--contact "{absolute_path(args.contact)}" '
        f'--raw-contact "{absolute_path(args.raw_contact)}" '
        f'--outdir "{output_dir}" '
        f'-t {args.threads}'
    )

    if args.genomad_db:
        command += f' --genomad_db "{absolute_path(args.genomad_db)}"'
    if args.checkv_db:
        command += f' --checkv_db "{absolute_path(args.checkv_db)}"'
    command += (
        f' --genomad-splits {args.genomad_splits}'
        f' --genomad-sensitivity {args.genomad_sensitivity}'
        f' --genomad-preset "{args.genomad_preset}"'
        f' --genomad-min-score {args.genomad_min_score}'
        f' --genomad-max-fdr {args.genomad_max_fdr}'
        f' --checkv-trimming "{args.checkv_trimming}"'
        f' --association-filter "{args.association_filter}"'
        f' --zscore-threshold {args.zscore_threshold}'
        f' --fixed-contact-threshold {args.fixed_contact_threshold}'
        f' --top-percent {args.top_percent}'
        f' --min-raw-contacts {args.min_raw_contacts}'
        f' --viral-quality-levels "{args.viral_quality_levels}"'
        f' --min-terminal-overlap {args.min_terminal_overlap}'
        f' --max-terminal-overlap {args.max_terminal_overlap}'
        f' --min-contact-strength {args.min_contact_strength}'
    )
    if args.genomad_cleanup:
        command += " --genomad-cleanup"
    else:
        command += " --no-genomad-cleanup"
    if args.genomad_restart:
        command += " --genomad-restart"
    if args.genomad_extra_args:
        command += f' --genomad-extra-args "{args.genomad_extra_args}"'
    if args.checkv_remove_tmp:
        command += " --checkv-remove-tmp"
    if args.checkv_restart:
        command += " --checkv-restart"
    if args.checkv_extra_args:
        command += f' --checkv-extra-args "{args.checkv_extra_args}"'
    if args.tmp_dir:
        command += f' --tmp-dir "{absolute_path(args.tmp_dir)}"'

    run_command(command)

def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "preprocessing" and sys.argv[2] == "--print-defaults":
        print(preprocessing_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "assembly" and sys.argv[2] == "--print-defaults":
        print(assembly_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "alignment" and sys.argv[2] == "--print-defaults":
        print(alignment_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "coverage" and sys.argv[2] == "--print-defaults":
        print(coverage_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "contact" and sys.argv[2] == "--print-defaults":
        print(contact_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "reassembly" and sys.argv[2] == "--print-defaults":
        print(reassembly_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "scaffolding" and sys.argv[2] == "--print-defaults":
        print(scaffolding_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "annotation" and sys.argv[2] == "--print-defaults":
        print(annotation_defaults_text())
        return
    if len(sys.argv) >= 3 and sys.argv[1] == "mge" and sys.argv[2] == "--print-defaults":
        print(mge_defaults_text())
        return

    parser = argparse.ArgumentParser(description="METAHICT Pipeline Wrapper")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Preprocessing subcommand 
    pre_parser = subparsers.add_parser("preprocessing", help="Run preprocessing of shotgun or Hi-C reads")
    pre_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    pre_parser.add_argument("-1", dest="r1", required=True, help="Forward shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)")
    pre_parser.add_argument("-2", dest="r2", required=True, help="Reverse shotgun or Hi-C FASTQ reads (.fastq or .fastq.gz)")
    pre_parser.add_argument("-o", "--output", required=True, help="Output directory for preprocessed reads")
    pre_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    pre_parser.add_argument("--dedup", action="store_true", help="Enable duplicate removal for Hi-C reads")
    pre_parser.add_argument("--prefix", help="Custom prefix for output files (default=derived from input filename)")
    pre_parser.add_argument("--minlen", type=int, default=50, help="Minimum read length after trimming (default=50)")
    pre_parser.add_argument("--trimq", type=int, default=10, help="Quality threshold for trimming (default=10)")
    pre_parser.add_argument("--qtrim", default="r", help="BBDuk quality trimming direction (default=r)")
    pre_parser.add_argument("--ftl", type=int, default=10, help="Trim bases from the left (default=10)")
    pre_parser.add_argument("--xmx", help="Java memory for BBDuk/BBMap tools (default=80%% of available memory)")
    pre_parser.add_argument("--ftm", type=int, default=5, help="BBDuk ftm value for modulo trimming (default=5)")
    pre_parser.add_argument("--ktrim", default="r", help="BBDuk adapter trimming direction (default=r)")
    pre_parser.add_argument("--k", type=int, default=23, help="BBDuk adapter k-mer size (default=23)")
    pre_parser.add_argument("--mink", type=int, default=11, help="BBDuk minimum adapter k-mer size (default=11)")
    pre_parser.add_argument("--hdist", type=int, default=1, help="BBDuk adapter k-mer Hamming distance (default=1)")
    pre_parser.add_argument("--adapter-ref", help="Adapter reference FASTA (default=<project_path>/external/bbmap/resources/adapters.fa)")
    pre_parser.add_argument("--skip-pre-qc-report", action="store_true", help="Skip FastQC report for input reads")
    pre_parser.add_argument("--skip-post-qc-report", action="store_true", help="Skip FastQC report for final reads")
    pre_parser.add_argument("--print-defaults", action="store_true", help="Print preprocessing defaults and exit")
    pre_parser.set_defaults(func=preprocessing)

    # Assembly subcommand
    asm_parser = subparsers.add_parser("assembly", help="Run assembly module to generate contigs")
    asm_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    asm_parser.add_argument("-1", dest="r1", required=True, help="Forward preprocessed reads (.fastq or .fastq.gz)")
    asm_parser.add_argument("-2", dest="r2", required=True, help="Reverse preprocessed reads (.fastq or .fastq.gz)")
    asm_parser.add_argument("-o", "--output", required=True, help="Output directory for assembled contigs")
    asm_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    asm_parser.add_argument("-m", "--memory", type=int, help="Memory in GB (default=80%% of available system memory)")
    asm_parser.add_argument("-l", "--min-len", type=int, default=1000, help="Minimum contig length (default=1000 bp)")
    asm_group = asm_parser.add_mutually_exclusive_group()
    asm_group.add_argument("--megahit", action="store_true", help="Assemble with MEGAHIT (default)")
    asm_group.add_argument("--metaspades", action="store_true", help="Assemble with metaSPAdes")
    asm_group.add_argument("--metaflye", action="store_true", help="Assemble with metaFlye")
    asm_parser.add_argument("--k-min", type=int, default=21, help="MEGAHIT minimum k-mer size (default=21)")
    asm_parser.add_argument("--k-max", type=int, default=141, help="MEGAHIT maximum k-mer size (default=141)")
    asm_parser.add_argument("--k-step", type=int, default=12, help="MEGAHIT k-mer step size (default=12)")
    asm_parser.add_argument("--merge-level", default="20,0.95", help="MEGAHIT merge level (default=20,0.95)")
    asm_parser.add_argument("--k-list", default="21,33,55", help="metaSPAdes k-mer sizes (default=21,33,55)")
    asm_parser.add_argument("--flye-method", default="--nano-raw", help="Flye read mode (default=--nano-raw)")
    asm_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    asm_parser.add_argument("--skip-quast", action="store_true", help="Skip QUAST QC report")
    asm_parser.add_argument("--keep-temp", action="store_true", help="Keep temporary assembler files")
    asm_parser.add_argument("--print-defaults", action="store_true", help="Print assembly defaults and exit")
    asm_parser.set_defaults(func=assembly)

    # Alignment subcommand
    aln_parser = subparsers.add_parser("alignment", help="Run alignment module for Hi-C read mapping")
    aln_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    aln_parser.add_argument("-r", "--reference", required=True, help="Assembled contigs file (.fasta)")
    aln_parser.add_argument("-1", dest="r1", required=True, help="Forward preprocessed Hi-C reads (.fastq or .fastq.gz)")
    aln_parser.add_argument("-2", dest="r2", required=True, help="Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)")
    aln_parser.add_argument("-o", "--output", required=True, help="Output directory for alignment results")
    aln_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    aln_parser.add_argument("--samtools-filter", default="-F 0x900", help="Filtering flag for samtools view (default='-F 0x900')")
    aln_parser.add_argument("--bwa-options", default="-5SP", help="BWA MEM options (default='-5SP')")
    aln_parser.add_argument("--mapq", type=int, default=30, help="Minimum MAPQ retained in BAM and used for metrics (default=30)")
    aln_parser.add_argument("--min-intra-dist", type=int, default=10000, help="Minimum same-contig distance for informative pairs (default=10000)")
    aln_parser.add_argument("--min-match-len", type=int, default=30, help="Minimum aligned nucleotide match length retained in BAM and used for metrics (default=30)")
    aln_parser.add_argument("--sort-memory", default="1G", help="Memory per samtools sort thread (default=1G)")
    aln_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    aln_parser.add_argument("--keep-sam", action="store_true", help="Keep intermediate map.sam")
    aln_parser.add_argument("--skip-metrics", action="store_true", help="Skip 3D ratio and informative-pair metrics")
    aln_parser.add_argument("--print-defaults", action="store_true", help="Print alignment defaults and exit")
    aln_parser.set_defaults(func=alignment)
    
    # Coverage subcommand
    cov_parser = subparsers.add_parser("coverage", help="Run coverage module for read mapping and coverage estimation")
    cov_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    cov_parser.add_argument("-1", dest="r1", required=True, help="Forward shotgun reads (.fastq or .fastq.gz)")
    cov_parser.add_argument("-2", dest="r2", required=True, help="Reverse shotgun reads (.fastq or .fastq.gz)")
    cov_parser.add_argument("-r", "--reference", required=True, help="Assembled contigs file (.fasta)")
    cov_parser.add_argument("-o", "--output", required=True, help="Output directory for coverage results")
    cov_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    cov_parser.add_argument("-m", "--memory", help="Java heap for BBMap, for example 64g or 64000m (default=80%% available memory)")
    cov_parser.add_argument("--percent-identity", type=float, default=97, help="Minimum end-to-end percent identity for coverage counting (default=97)")
    cov_parser.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality for coverage counting (default=0)")
    cov_parser.add_argument("--weight-mapq", type=float, default=0.0, help="Weight per-base depth by mapping quality (default=0.0; disabled)")
    cov_parser.add_argument("--include-edge-bases", action="store_true", help="Include read-length edge bases in depth and variance calculation")
    cov_parser.add_argument("--max-edge-bases", type=int, default=75, help="Maximum edge length excluded when edge bases are not included (default=75)")
    cov_parser.add_argument("--min-contig-length", type=int, default=0, help="Minimum contig length included by JGI coverage summarization (default=0; disabled)")
    cov_parser.add_argument("--min-contig-depth", type=float, default=0, help="Minimum contig depth used by JGI coverage summarization (default=0; disabled)")
    cov_parser.add_argument("--bbmap-extra-args", default="", help="Additional options passed directly to bbmap.sh (default=empty)")
    cov_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    cov_parser.add_argument("--keep-sam", action="store_true", help="Keep intermediate SG_map.sam")
    cov_parser.add_argument("--keep-temp", action="store_true", help="Keep temporary working/index files")
    cov_parser.add_argument("--print-defaults", action="store_true", help="Print coverage defaults and exit")
    cov_parser.set_defaults(func=coverage)
    
    # Contact subcommand
    con_parser = subparsers.add_parser("contact", help="Run contact module to generate contact matrices")
    con_parser.add_argument("method", choices=["raw", "normcc", "hiczin", "bin3c", "metator"], help="Normalization method")
    con_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    con_parser.add_argument("--bam", required=True, help="Hi-C read alignment file (.bam)")
    con_parser.add_argument("--fasta", required=True, help="Assembled contigs file (.fasta)")
    con_parser.add_argument("--out", dest="output", required=True, help="Output directory for contact maps")
    con_parser.add_argument("--enzyme", required=True, help="Restriction enzymes used in Hi-C library (e.g., Sau3AI,MluCI)")
    con_parser.add_argument("--metacc-min-signal", type=int, default=1, help="Minimum off-diagonal contact signal for retained contigs (default=1)")
    con_parser.add_argument("--metacc-min-len", type=int, default=1000, help="Minimum contig length for raw contact generation (default=1000)")
    con_parser.add_argument("--metacc-min-mapq", type=int, default=30, help="Minimum MAPQ for raw contact generation (default=30)")
    con_parser.add_argument("--metacc-min-match", type=int, default=30, help="Minimum terminal aligned match length for raw contact generation (default=30)")
    con_parser.add_argument("--spurious-contact-percent", type=float, default=5, help="Percentile cutoff for spurious-contact denoising (default=5)")
    con_parser.add_argument("--coverage-file", help="Coverage file from Module 4; required for hiczin and metator")
    con_parser.add_argument("--epsilon", type=float, default=1, help="Epsilon used by hiczin, bin3c, and metator normalization (default=1)")
    con_parser.add_argument("--max-iter", type=int, default=1000, help="Maximum iterations for bin3c matrix balancing (default=1000)")
    con_parser.add_argument("--tol", type=float, default=1e-6, help="Convergence tolerance for bin3c matrix balancing (default=1e-6)")
    con_parser.add_argument("--print-defaults", action="store_true", help="Print contact defaults and exit")
    con_parser.set_defaults(func=contact)
    
    # Binning subcommand
    bin_parser = subparsers.add_parser("binning", help="Run binning module to generate MAGs and refined bins")
    bin_parser.add_argument("--fasta", help="Assembled contigs file (.fa or .fasta)")
    bin_parser.add_argument("--bam", help="Hi-C reads aligned to the contigs (.bam)")
    bin_parser.add_argument("--output", help="Output directory for binning results")
    bin_parser.add_argument("--project_path", help="Path to the METAHICT project directory")
    bin_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    bin_parser.add_argument("--checkm_db", help="Path to CheckM database (if not using default environment variable)")
    bin_parser.add_argument("--enzyme", default="Sau3AI,MluCI", help="Comma-separated restriction enzymes (default=Sau3AI,MluCI)")
    bin_parser.add_argument("--metacc-min-len", type=int, default=1000, help="Minimum contig length for MetaCC contact generation (default=1000)")
    bin_parser.add_argument("--metacc-min-signal", type=int, default=2, help="Minimum Hi-C signal for MetaCC contact generation (default=2)")
    bin_parser.add_argument("--metacc-min-mapq", type=int, default=30, help="Minimum MAPQ for MetaCC contact generation (default=30)")
    bin_parser.add_argument("--metacc-min-match", type=int, default=30, help="Minimum aligned match length for MetaCC contact generation (default=30)")
    bin_parser.add_argument("--metacc-min-binsize", type=int, default=150000, help="Minimum MetaCC output bin size (default=150000)")
    bin_parser.add_argument("--normcc-thres", type=float, default=0.05, help="Fraction of NormCC-normalized contacts discarded (default=0.05)")
    bin_parser.add_argument("--bin3c-min-len", type=int, default=1000, help="Minimum contig length for bin3C-compatible contact generation (default=1000)")
    bin_parser.add_argument("--bin3c-min-signal", type=int, default=5, help="Minimum Hi-C signal for bin3C-compatible contact generation (default=5)")
    bin_parser.add_argument("--bin3c-min-mapq", type=int, default=60, help="Minimum MAPQ for bin3C-compatible contact generation (default=60)")
    bin_parser.add_argument("--bin3c-min-match", type=int, default=10, help="Minimum aligned match length for bin3C-compatible contact generation (default=10)")
    bin_parser.add_argument("--bin3c-min-extent", type=int, default=50000, help="Minimum bin3C cluster extent (default=50000)")
    bin_parser.add_argument("--min-completeness", type=float, default=50, help="Final METAHICT minimum completeness (default=50)")
    bin_parser.add_argument("--max-contamination", type=float, default=10, help="Final METAHICT maximum contamination (default=10)")
    bin_parser.add_argument("--contamination-penalty", type=float, default=5, help="Penalty used in completeness - penalty * contamination (default=5)")
    bin_parser.add_argument("--min-input-bin-size", type=int, default=50000, help="Minimum input bin FASTA file size before refinement (default=50000 bytes)")
    bin_parser.add_argument("--max-input-bin-size", type=int, default=20000000, help="Maximum input bin FASTA file size before refinement (default=20000000 bytes)")
    bin_parser.add_argument("--binning-refiner-min-size", type=int, default=524288, help="Minimum refined bin size for Binning_refiner (default=524288 bp)")
    bin_parser.add_argument("--tmp-dir", default=None, help="Temporary directory root for CheckM2 working files (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    bin_parser.add_argument("--keep-temp", action="store_true", help="Keep successful CheckM2 temporary directories")
    bin_parser.add_argument("--no-fasta", action="store_true", help="Do not write bin3C cluster FASTA files")
    bin_parser.add_argument("--no-report", action="store_true", help="Do not write bin3C cluster report")
    bin_parser.add_argument("--no-spades", action="store_true", help="Input assembly was not produced by SPAdes")
    bin_parser.add_argument("--only-large", action="store_true", help="Only write bin3C FASTA clusters larger than --bin3c-min-extent")
    bin_parser.add_argument("--skip-checkm2", action="store_true", help="Skip CheckM2 during final bin refinement")
    bin_parser.add_argument("--skip-refinement", action="store_true", help="Skip Binning_refiner combinations")
    bin_parser.add_argument("--skip-consolidation", action="store_true", help="Skip final consolidation across bin sets")
    bin_parser.add_argument("--keep-ambiguous", action="store_true", help="Keep ambiguous contigs in all bins")
    bin_parser.add_argument("--remove-ambiguous", action="store_true", help="Remove ambiguous contigs from all bins")
    bin_parser.add_argument("--seed", type=int, help="Random seed (default=empty)")
    bin_parser.add_argument("--print-defaults", action="store_true", help="Print binning defaults and exit")
    bin_parser.set_defaults(func=binning)
    
    # Reassembly subcommand
    reas_parser = subparsers.add_parser("reassembly", help="Run reassembly module to reconstruct bins and unmapped reads")
    reas_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    reas_parser.add_argument("--bin", required=True, help="Directory containing input bins")
    reas_parser.add_argument("--assembly", required=True, help="Original assembly FASTA file (.fa or .fasta)")
    reas_parser.add_argument("--hic1", required=True, help="Forward preprocessed Hi-C reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--hic2", required=True, help="Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--sg1", required=True, help="Forward preprocessed shotgun reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--sg2", required=True, help="Reverse preprocessed shotgun reads (.fastq or .fastq.gz)")
    reas_parser.add_argument("--bam", required=True, help="Hi-C read alignments to the assembly (.bam)")
    reas_parser.add_argument("--outdir", required=True, help="Output directory for reassembly results")
    reas_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    reas_parser.add_argument("-m", "--memory", type=int, help="Memory in GB (default=80%% of available system memory)")
    reas_parser.add_argument("--cutoff-quantile", type=float, default=0.95, help="Quantile of the short-insert component used for selection (default=0.95)")
    reas_parser.add_argument("--top-k", type=int, default=100, help="Number of longest contigs used for EM fitting (default=100)")
    reas_parser.add_argument("--min-mapq", type=int, default=30, help="Minimum MAPQ for pair-level insert-size extraction (default=30)")
    reas_parser.add_argument("--min-match-len", type=int, default=30, help="Minimum aligned match length for pair-level extraction (default=30)")
    reas_parser.add_argument("--exclude-duplicates", action="store_true", help="Exclude duplicate-marked alignments")
    reas_parser.add_argument("--write-nonselected-hic", action="store_true", help="Also write non-selected Hi-C FASTQ files")
    reas_parser.add_argument("--min-contig-len", type=int, default=500, help="Minimum reassembled contig length retained (default=500)")
    reas_parser.add_argument("--strict-cut-off", type=int, default=2, help="Strict read-recruitment SNP cutoff (default=2)")
    reas_parser.add_argument("--permissive-cut-off", type=int, default=5, help="Permissive read-recruitment SNP cutoff (default=5)")
    reas_parser.add_argument("--contamination-penalty", type=float, default=5, help="Penalty used in completeness - penalty * contamination (default=5)")
    reas_parser.add_argument("--skip-checkm2", action="store_true", help="Skip CheckM2 quality evaluation")
    reas_parser.add_argument("--checkm2_db", help="Path to CheckM2 database (if not using default environment variable)")
    reas_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    reas_parser.add_argument("--spades-mode", choices=["careful", "none"], default="careful", help="SPAdes mode (default=careful)")
    reas_parser.add_argument("--spades-phred-offset", help="SPAdes phred offset (default=empty)")
    reas_parser.add_argument("--spades-extra-args", default="", help="Additional options passed directly to SPAdes (default=empty)")
    reas_parser.add_argument("--skip-residual-assembly", action="store_true", help="Do not assemble residual unmapped reads")
    reas_parser.add_argument("--keep-temp", action="store_true", help="Keep successful SPAdes and CheckM2 temporary directories")
    reas_parser.add_argument("--print-defaults", action="store_true", help="Print reassembly defaults and exit")
    reas_parser.set_defaults(func=reassembly)

    # Scaffolding subcommand
    scaf_parser = subparsers.add_parser("scaffolding", help="Run scaffolding module to generate scaffolded genomes and heatmaps")
    scaf_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    scaf_parser.add_argument("--fasta", required=True, help="Input bin FASTA file for scaffolding (.fa or .fasta)")
    scaf_parser.add_argument("--bam", help="Optional Hi-C read alignments to the assembly (.bam)")
    scaf_parser.add_argument("--enzyme", required=True, help="Restriction enzymes used in the Hi-C library (e.g., Sau3AI,MluCI)")
    scaf_parser.add_argument("--outdir", required=True, help="Output directory for scaffolding results")
    scaf_parser.add_argument("--hic1", required=True, help="Forward preprocessed Hi-C reads (.fastq or .fastq.gz)")
    scaf_parser.add_argument("--hic2", required=True, help="Reverse preprocessed Hi-C reads (.fastq or .fastq.gz)")
    scaf_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    scaf_parser.add_argument("-m", "--memory", type=str, default=None, help="Memory limit (default=80%% of available RAM)")
    scaf_parser.add_argument("-r", "--resolution", type=int, default=10000, help="Segment length for visualization (default=10000 bp)")
    scaf_parser.add_argument("--min-contig-len", type=int, default=5000, help="Minimum contig length retained for scaffolding (default=5000)")
    scaf_parser.add_argument("--bwa-options", default="-5SP", help="BWA MEM options for internal alignments (default='-5SP')")
    scaf_parser.add_argument("--samtools-filter", default="-F 0x900", help="samtools view filter for internal alignments (default='-F 0x900')")
    scaf_parser.add_argument("--sort-memory", default="1G", help="Memory per samtools sort thread for internal alignments (default=1G)")
    scaf_parser.add_argument("--metacc-min-mapq", type=int, default=30, help="Minimum MAPQ for MetaCC contact generation (default=30)")
    scaf_parser.add_argument("--metacc-min-len", type=int, default=1000, help="Minimum contig length for MetaCC contact generation (default=1000)")
    scaf_parser.add_argument("--metacc-min-match", type=int, default=30, help="Minimum aligned match length for MetaCC contact generation (default=30)")
    scaf_parser.add_argument("--metacc-min-signal", type=int, default=2, help="Minimum signal for MetaCC contact generation (default=2)")
    scaf_parser.add_argument("--bin3c-min-mapq", type=int, default=60, help="Minimum MAPQ for bin3C-compatible contact generation (default=60)")
    scaf_parser.add_argument("--bin3c-min-len", type=int, default=1000, help="Minimum contig length for bin3C-compatible contact generation (default=1000)")
    scaf_parser.add_argument("--bin3c-min-match", type=int, default=10, help="Minimum aligned match length for bin3C-compatible contact generation (default=10)")
    scaf_parser.add_argument("--bin3c-min-signal", type=int, default=5, help="Minimum signal for bin3C-compatible contact generation (default=5)")
    scaf_parser.add_argument("--yahs-resolutions", default="", help="YaHS scaffolding resolution list; empty uses YaHS automatic selection (default=empty)")
    scaf_parser.add_argument("--yahs-min-mapq", type=int, default=10, help="YaHS minimum mapping quality (default=10)")
    scaf_parser.add_argument("--yahs-min-contig-len", type=int, default=0, help="YaHS minimum contig length to scaffold (default=0)")
    scaf_parser.add_argument("--yahs-rounds", type=int, default=1, help="YaHS rounds at each resolution level (default=1)")
    scaf_parser.add_argument("--yahs-no-contig-ec", action="store_true", help="Disable YaHS contig error correction")
    scaf_parser.add_argument("--yahs-no-scaffold-ec", action="store_true", help="Disable YaHS scaffold error correction")
    scaf_parser.add_argument("--yahs-no-mem-check", action="store_true", help="Disable YaHS runtime memory check")
    scaf_parser.add_argument("--yahs-extra-args", default="", help="Additional options passed directly to YaHS (default=empty)")
    scaf_parser.add_argument("--normcc-thres", type=float, default=0.05, help="NormCC denoising threshold (default=0.05)")
    scaf_parser.add_argument("--heatmap-max-image", type=int, default=5000, help="Maximum heatmap image dimension before downsampling (default=5000)")
    scaf_parser.add_argument("--skip-checkm2", action="store_true", help="Skip CheckM2 quality evaluation")
    scaf_parser.add_argument("--checkm2_db", help="Path to CheckM2 database (if not using default environment variable)")
    scaf_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    scaf_parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files for debugging")
    scaf_parser.add_argument("--print-defaults", action="store_true", help="Print scaffolding defaults and exit")
    scaf_parser.set_defaults(func=scaffolding)

    # Annotation subcommand
    anno_parser = subparsers.add_parser("annotation", help="Run annotation module using GTDB-Tk for taxonomic classification")
    anno_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    anno_parser.add_argument("--bin", required=True, help="Directory containing input bins")
    anno_parser.add_argument("--outdir", required=True, help="Output directory for annotation results")
    anno_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    anno_parser.add_argument("--pplacer-cpus", type=int, help="Number of CPUs for pplacer (default=same as --threads)")
    anno_parser.add_argument("--gtdbtk_db", "--gtdbtk-db", dest="gtdbtk_db", help="Path to GTDB-Tk database (default=<project_path>/databases/gtdbtk_db/release220)")
    anno_parser.add_argument("--extension", default="fa", help="Genome file extension processed by GTDB-Tk (default=fa)")
    anno_parser.add_argument("--prefix", default="gtdbtk", help="Prefix for GTDB-Tk output files (default=gtdbtk)")
    anno_parser.add_argument("--skip-ani-screen", dest="skip_ani_screen", action="store_true", help="Skip GTDB-Tk ANI screening (default=true)")
    anno_parser.add_argument("--no-skip-ani-screen", dest="skip_ani_screen", action="store_false", help="Use GTDB-Tk ANI/Mash screening; requires --mash-db")
    anno_parser.set_defaults(skip_ani_screen=True)
    anno_parser.add_argument("--mash-db", help="Mash database path used when ANI screening is enabled")
    anno_parser.add_argument("--min-perc-aa", type=float, default=10, help="Minimum percentage of amino acids in the MSA (default=10)")
    anno_parser.add_argument("--min-af", type=float, default=0.5, help="Minimum alignment fraction for species assignment (default=0.5)")
    anno_parser.add_argument("--full-tree", action="store_true", help="Use GTDB-Tk full bacterial tree (default=false)")
    anno_parser.add_argument("--scratch-dir", help="Scratch directory for pplacer disk-backed memory reduction (default=empty)")
    anno_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    anno_parser.add_argument("--force", action="store_true", help="Continue if GTDB-Tk errors on a single genome")
    anno_parser.add_argument("--keep-intermediates", action="store_true", help="Keep GTDB-Tk intermediate files")
    anno_parser.add_argument("--debug", action="store_true", help="Keep GTDB-Tk debug intermediates")
    anno_parser.add_argument("--write-single-copy-genes", action="store_true", help="Output unaligned single-copy marker genes")
    anno_parser.add_argument("--gtdbtk-extra-args", default="", help='Additional native options passed to gtdbtk classify_wf, for example "--genes" or "--no_mash --mash_k 21" (default=empty)')
    anno_parser.add_argument("--print-defaults", action="store_true", help="Print annotation defaults and exit")
    anno_parser.set_defaults(func=annotation)

    # MGE subcommand
    mge_parser = subparsers.add_parser("mge", help="Run MGE module for viral/plasmid detection and candidate MGE-MAG association analysis")
    mge_parser.add_argument("-p", "--project_path", required=True, help="Path to the METAHICT project directory")
    mge_parser.add_argument("--combined", required=True, help="Combined contigs FASTA file (includes binned and unmapped contigs)")
    mge_parser.add_argument("--contact", required=True, help="Normalized contact matrix (.npz)")
    mge_parser.add_argument("--raw-contact", dest="raw_contact", required=True, help="Raw contact matrix (.npz) used for raw Hi-C support")
    mge_parser.add_argument("--outdir", required=True, help="Output directory for MGE analysis results")
    mge_parser.add_argument("-t", "--threads", type=int, default=80, help="Number of CPU threads (default=80)")
    mge_parser.add_argument("--genomad_db", "--genomad-db", dest="genomad_db", help="Path to geNomad database (default=<project_path>/databases/genomad_db/genomad_db)")
    mge_parser.add_argument("--checkv_db", "--checkv-db", dest="checkv_db", help="Path to CheckV database (default=<project_path>/databases/checkv_db/checkv-db-v1.5)")
    mge_parser.add_argument("--genomad-splits", type=int, default=8, help="geNomad MMseqs2 split count (default=8)")
    mge_parser.add_argument("--genomad-sensitivity", type=float, default=4.2, help="geNomad MMseqs2 sensitivity (default=4.2)")
    mge_parser.add_argument("--genomad-cleanup", dest="genomad_cleanup", action="store_true", help="Delete geNomad intermediate files (default=true)")
    mge_parser.add_argument("--no-genomad-cleanup", dest="genomad_cleanup", action="store_false", help="Keep geNomad intermediate files")
    mge_parser.set_defaults(genomad_cleanup=True)
    mge_parser.add_argument("--genomad-restart", action="store_true", help="Overwrite existing geNomad intermediate files")
    mge_parser.add_argument("--genomad-preset", choices=["default", "conservative", "relaxed"], default="default", help="geNomad filtering preset (default=default)")
    mge_parser.add_argument("--genomad-min-score", type=float, default=0.7, help="geNomad minimum virus/plasmid score when preset is default (default=0.7)")
    mge_parser.add_argument("--genomad-max-fdr", type=float, default=0.1, help="geNomad maximum FDR when preset is default (default=0.1)")
    mge_parser.add_argument("--genomad-extra-args", default="", help="Additional native options passed to genomad end-to-end (default=empty)")
    mge_parser.add_argument("--checkv-remove-tmp", action="store_true", help="Delete CheckV intermediate files")
    mge_parser.add_argument("--checkv-restart", action="store_true", help="Overwrite existing CheckV intermediate files")
    mge_parser.add_argument("--checkv-trimming", choices=["default", "conservative", "aggressive"], default="default", help="CheckV trimming mode (default=default)")
    mge_parser.add_argument("--checkv-extra-args", default="", help="Additional native options passed to checkv end_to_end (default=empty)")
    mge_parser.add_argument("--association-filter", choices=["zscore", "fixed", "percentage", "raw-support-only"], default="zscore", help="Association filter method (default=zscore)")
    mge_parser.add_argument("--zscore-threshold", type=float, default=0.5, help="Minimum Z-score for filtered viral MGE-MAG associations (default=0.5)")
    mge_parser.add_argument("--fixed-contact-threshold", type=float, default=0, help="Minimum contact strength for fixed association filtering (default=0)")
    mge_parser.add_argument("--top-percent", type=float, default=50, help="Top percent of association contact strengths retained by percentage filtering (default=50)")
    mge_parser.add_argument("--min-raw-contacts", type=float, default=2, help="Minimum raw Hi-C read-pair support for candidate associations (default=2)")
    mge_parser.add_argument("--viral-quality-levels", default="Complete,High-quality,Medium-quality,Low-quality", help="Comma-separated CheckV quality levels used for candidate viral MGE-MAG association")
    mge_parser.add_argument("--min-terminal-overlap", type=int, default=50, help="Minimum terminal overlap length for circularity evidence (default=50)")
    mge_parser.add_argument("--max-terminal-overlap", type=int, default=1000, help="Maximum terminal overlap length for circularity evidence (default=1000)")
    mge_parser.add_argument("--min-contact-strength", type=float, default=0, help="Minimum positive contact strength counted for viral MGE-MAG associations (default=0)")
    mge_parser.add_argument("--tmp-dir", help="Temporary directory root (default=METAHICT_TMP_ROOT, TMPDIR, or /tmp)")
    mge_parser.add_argument("--print-defaults", action="store_true", help="Print MGE defaults and exit")
    mge_parser.set_defaults(func=mge)

    args = parser.parse_args()
    if args.command == "preprocessing":
        preprocessing(args)
        
    elif args.command == "assembly":
        assembly(args)
        
    elif args.command == "alignment":
        alignment(args)

    elif args.command == "coverage":
        coverage(args)
        
    elif args.command == "contact":
        contact(args)
        
    elif args.command == "binning":
        binning(args)
        
    elif args.command == "reassembly":
        reassembly(args)

    elif args.command == "scaffolding":
        scaffolding(args)

    elif args.command == "annotation":
        annotation(args)

    elif args.command == "mge":
        mge(args)

if __name__ == "__main__":
    main()
