/*
 * Native DSL2 stage processes.  Each process is bound to one METAHICT module
 * and declares the material files it consumes and produces.  The module
 * implementation remains the scientific reference code; Nextflow owns the
 * dependency graph, task resources, environments, and published outputs.
 */

process CONDA_BUNDLE {
    tag 'conda-environment-bundle'
    cpus 1
    memory '1 GB'
    time '30m'

    input:
    val project_path

    output:
    path 'conda_bundle.ok', emit: ready

    script:
    """
    set -euo pipefail
    bash "${project_path}/installation/verify_env_bundle.sh" --project "${project_path}"
    touch conda_bundle.ok
    """
}

process PREPROCESSING {
    tag "${sample}:${library}"
    cpus 8
    memory '16 GB'
    time '12h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.preprocessing }
    publishDir { "${params.out_root}/${sample}/1_preprocessing/${library}" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), val(library), path(read1), path(read2)

    output:
    tuple val(sample), val(row), val(library), path('preprocessing'), emit: results

    script:
    def extra = row["preprocessing_${library}_extra_args"] ?: ''
    def dedup = library == 'hic' ? '--dedup' : ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    mkdir -p preprocessing
    python "${params.project_path}/metahict.py" preprocessing \\
      -p "${params.project_path}" -1 "${read1}" -2 "${read2}" \\
      -o preprocessing -t "${task.cpus}" ${dedup} ${extra}
    test -n "\$(find -L preprocessing -name 'final_*_1.fastq.gz' -print -quit)"
    test -n "\$(find -L preprocessing -name 'final_*_2.fastq.gz' -print -quit)"
    """

    stub:
    """
    mkdir -p preprocessing
    touch preprocessing/final_${library}_1.fastq.gz preprocessing/final_${library}_2.fastq.gz
    """
}

process ASSEMBLY {
    tag "${sample}:assembly"
    cpus 32
    memory '128 GB'
    time '48h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.assembly }
    publishDir { "${params.out_root}/${sample}/2_assembly" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(sg_dir)

    output:
    tuple val(sample), val(row), path('assembly'), emit: results

    script:
    def extra = row.assembly_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    sg1=\$(find -L "${sg_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    sg2=\$(find -L "${sg_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -n "\${sg1}" && test -n "\${sg2}"
    mkdir -p assembly
    python "${params.project_path}/metahict.py" assembly \\
      -p "${params.project_path}" -1 "\${sg1}" -2 "\${sg2}" \\
      -o assembly -t "${task.cpus}" --megahit ${extra}
    test -s assembly/final_assembly.fasta
    """

    stub:
    """
    mkdir -p assembly
    touch assembly/final_assembly.fasta
    """
}

process ALIGNMENT {
    tag "${sample}:alignment"
    cpus 24
    memory '64 GB'
    time '24h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.alignment }
    publishDir { "${params.out_root}/${sample}/3_alignment" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(assembly_dir), path(hic_dir)

    output:
    tuple val(sample), val(row), path('alignment'), emit: results

    script:
    def extra = row.alignment_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -n "\${hic1}" && test -n "\${hic2}"
    mkdir -p alignment
    python "${params.project_path}/metahict.py" alignment \\
      -p "${params.project_path}" -r "${assembly_dir}/final_assembly.fasta" \\
      -1 "\${hic1}" -2 "\${hic2}" -o alignment -t "${task.cpus}" ${extra}
    test -s alignment/sorted_map.bam
    """

    stub:
    """
    mkdir -p alignment
    touch alignment/sorted_map.bam
    """
}

process COVERAGE {
    tag "${sample}:coverage"
    cpus 24
    memory '64 GB'
    time '24h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.coverage }
    publishDir { "${params.out_root}/${sample}/4_coverage" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(assembly_dir), path(sg_dir)

    output:
    tuple val(sample), val(row), path('coverage'), emit: results

    script:
    def extra = row.coverage_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    sg1=\$(find -L "${sg_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    sg2=\$(find -L "${sg_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -n "\${sg1}" && test -n "\${sg2}"
    mkdir -p coverage
    python "${params.project_path}/metahict.py" coverage \\
      -p "${params.project_path}" -r "${assembly_dir}/final_assembly.fasta" \\
      -1 "\${sg1}" -2 "\${sg2}" -o coverage -t "${task.cpus}" ${extra}
    test -s coverage/coverage.txt
    """

    stub:
    """
    mkdir -p coverage
    touch coverage/coverage.txt coverage/pair.txt
    """
}

process CONTACT {
    tag "${sample}:contact"
    cpus 16
    memory '48 GB'
    time '24h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.contact }
    publishDir { "${params.out_root}/${sample}/5_contact" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(assembly_dir), path(alignment_dir), path(coverage_dir)

    output:
    tuple val(sample), val(row), path('contact'), emit: results

    script:
    def extra = row.contact_extra_args ?: ''
    def enzyme = row.enzyme ?: 'Sau3AI,MluCI'
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    mkdir -p contact
    python "${params.project_path}/metahict.py" contact normcc \\
      -p "${params.project_path}" --bam "${alignment_dir}/sorted_map.bam" \\
      --fasta "${assembly_dir}/final_assembly.fasta" --out contact --enzyme "${enzyme}" ${extra}
    test -s contact/Raw_contact_matrix.npz
    """

    stub:
    """
    mkdir -p contact
    touch contact/Raw_contact_matrix.npz contact/denoised_contact_matrix_normcc.npz
    """
}

process BINNING {
    tag "${sample}:binning"
    cpus 32
    memory '128 GB'
    time '72h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.binning }
    publishDir { "${params.out_root}/${sample}/6_binning" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(assembly_dir), path(alignment_dir)

    output:
    tuple val(sample), val(row), path('binning'), emit: results

    script:
    def extra = row.binning_extra_args ?: ''
    def enzyme = row.enzyme ?: 'Sau3AI,MluCI'
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    export CHECKM2DB="${params.checkm2_db}"
    test -d "${params.checkm_db}"
    test -s "${params.checkm2_db}"
    mkdir -p binning
    python "${params.project_path}/metahict.py" binning \\
      --project_path "${params.project_path}/modules" \\
      --fasta "${assembly_dir}/final_assembly.fasta" --bam "${alignment_dir}/sorted_map.bam" \\
      --output binning --enzyme "${enzyme}" -t "${task.cpus}" --no-spades \\
      --checkm_db "${params.checkm_db}" ${extra}
    test -d binning/metahict/metahict_50_10_bins
    """

    stub:
    """
    mkdir -p binning/metahict/metahict_50_10_bins
    touch binning/metahict/metahict_50_10_bins/bin.1.fa
    """
}

process REASSEMBLY {
    tag "${sample}:reassembly"
    cpus 32
    memory '128 GB'
    time '72h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.reassembly }
    publishDir { "${params.out_root}/${sample}/7_reassembly" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(assembly_dir), path(alignment_dir), path(binning_dir), path(sg_dir, stageAs: 'sg_preprocessing'), path(hic_dir, stageAs: 'hic_preprocessing')

    output:
    tuple val(sample), val(row), path('reassembly'), emit: results

    script:
    def extra = row.reassembly_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    export CHECKM2DB="${params.checkm2_db}"
    test -s "${params.checkm2_db}"
    sg1=\$(find -L "${sg_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    sg2=\$(find -L "${sg_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    bins="${binning_dir}/metahict/metahict_50_10_bins"
    test -n "\${sg1}" && test -n "\${sg2}" && test -n "\${hic1}" && test -n "\${hic2}"
    mkdir -p reassembly
    python "${params.project_path}/metahict.py" reassembly \\
      -p "${params.project_path}/modules" --bin "\${bins}" \\
      --assembly "${assembly_dir}/final_assembly.fasta" --hic1 "\${hic1}" --hic2 "\${hic2}" \\
      --sg1 "\${sg1}" --sg2 "\${sg2}" --bam "${alignment_dir}/sorted_map.bam" \\
      --outdir reassembly -t "${task.cpus}" --checkm2_db "${params.checkm2_db}" ${extra}
    test -d reassembly/reassembled_bins
    """

    stub:
    """
    mkdir -p reassembly/reassembled_bins
    touch reassembly/reassembled_bins/bin1.fa reassembly/combined_contigs.fa
    """
}

process ANNOTATION {
    tag "${sample}:annotation"
    cpus 16
    memory '128 GB'
    time '48h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.annotation }
    publishDir { "${params.out_root}/${sample}/9_annotation" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(reassembly_dir)

    output:
    tuple val(sample), val(row), path('annotation'), emit: results

    script:
    def extra = row.annotation_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/gtdbtk-2.4.0/bin:\$PATH"
    test -d "${params.gtdbtk_db}"
    mkdir -p annotation
    python "${params.project_path}/metahict.py" annotation \\
      -p "${params.project_path}" --bin "${reassembly_dir}/reassembled_bins" \\
      --outdir annotation -t "${task.cpus}" --gtdbtk_db "${params.gtdbtk_db}" ${extra}
    test -s annotation/classify/gtdbtk.bac120.summary.tsv
    """

    stub:
    """
    mkdir -p annotation/classify
    touch annotation/classify/gtdbtk.bac120.summary.tsv
    """
}

process SCAFFOLDING {
    tag "${sample}:scaffolding"
    cpus 32
    memory '128 GB'
    time '72h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.scaffolding }
    publishDir { "${params.out_root}/${sample}/8_scaffolding" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(reassembly_dir), path(alignment_dir), path(hic_dir)

    output:
    tuple val(sample), val(row), path('scaffolding'), emit: results

    script:
    def extra = row.scaffolding_extra_args ?: ''
    def enzyme = row.enzyme ?: 'Sau3AI,MluCI'
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    export CHECKM2DB="${params.checkm2_db}"
    test -s "${params.checkm2_db}"
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    bin_fasta=\$(find -L "${reassembly_dir}/reassembled_bins" -type f \\( -name '*.fa' -o -name '*.fasta' \\) -print -quit)
    test -n "\${hic1}" && test -n "\${hic2}" && test -n "\${bin_fasta}"
    mkdir -p scaffolding
    python "${params.project_path}/metahict.py" scaffolding \\
      -p "${params.project_path}" --fasta "\${bin_fasta}" --enzyme "${enzyme}" \\
      --outdir scaffolding --hic1 "\${hic1}" --hic2 "\${hic2}" \\
      --bam "${alignment_dir}/sorted_map.bam" -t "${task.cpus}" \\
      --checkm2_db "${params.checkm2_db}" ${extra}
    test -d scaffolding
    """

    stub:
    """
    mkdir -p scaffolding
    touch scaffolding/scaffolds.fa
    """
}

process MGE_ALIGNMENT {
    tag "${sample}:mge-alignment"
    cpus 24
    memory '64 GB'
    time '24h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.alignment }
    publishDir { "${params.out_root}/${sample}/10_MGE/mge_alignment" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(reassembly_dir), path(hic_dir)

    output:
    tuple val(sample), val(row), path('mge_alignment'), emit: results

    script:
    def extra = row.mge_alignment_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -s "${reassembly_dir}/combined_contigs.fa"
    mkdir -p mge_alignment
    python "${params.project_path}/metahict.py" alignment \\
      -p "${params.project_path}" -r "${reassembly_dir}/combined_contigs.fa" \\
      -1 "\${hic1}" -2 "\${hic2}" -o mge_alignment -t "${task.cpus}" ${extra}
    test -s mge_alignment/sorted_map.bam
    """

    stub:
    """
    mkdir -p mge_alignment
    touch mge_alignment/sorted_map.bam
    """
}

process MGE_CONTACT {
    tag "${sample}:mge-contact"
    cpus 16
    memory '48 GB'
    time '24h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.contact }
    publishDir { "${params.out_root}/${sample}/10_MGE/mge_contact" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(reassembly_dir), path(mge_alignment_dir)

    output:
    tuple val(sample), val(row), path('mge_contact'), emit: results

    script:
    def extra = row.mge_contact_extra_args ?: ''
    def enzyme = row.enzyme ?: 'Sau3AI,MluCI'
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    mkdir -p mge_contact
    python "${params.project_path}/metahict.py" contact normcc \\
      -p "${params.project_path}" --bam "${mge_alignment_dir}/sorted_map.bam" \\
      --fasta "${reassembly_dir}/combined_contigs.fa" --out mge_contact --enzyme "${enzyme}" ${extra}
    test -s mge_contact/Raw_contact_matrix.npz
    test -s mge_contact/denoised_contact_matrix_normcc.npz
    """

    stub:
    """
    mkdir -p mge_contact
    touch mge_contact/Raw_contact_matrix.npz mge_contact/denoised_contact_matrix_normcc.npz
    """
}

process MGE {
    tag "${sample}:mge"
    cpus 16
    memory '64 GB'
    time '72h'
    conda params.metahict_conda
    container { params.container_image_override ?: params.container_images.mge }
    publishDir { "${params.out_root}/${sample}/10_MGE" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(row), path(reassembly_dir), path(mge_contact_dir)

    output:
    tuple val(sample), val(row), path('mge'), emit: results

    script:
    def extra = row.mge_extra_args ?: ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/genomad/bin:${params.conda_envs_path}/checkv_env/bin:\$PATH"
    test -d "${params.genomad_db}"
    test -d "${params.checkv_db}"
    mkdir -p mge
    python "${params.project_path}/metahict.py" mge \\
      -p "${params.project_path}" --combined "${reassembly_dir}/combined_contigs.fa" \\
      --contact "${mge_contact_dir}/denoised_contact_matrix_normcc.npz" \\
      --raw-contact "${mge_contact_dir}/Raw_contact_matrix.npz" \\
      --outdir mge -t "${task.cpus}" --genomad_db "${params.genomad_db}" \\
      --checkv_db "${params.checkv_db}" ${extra}
    test -d mge
    """

    stub:
    """
    mkdir -p mge
    touch mge/candidate_mge_mag_associations_zscore_filtered.tsv
    """
}
