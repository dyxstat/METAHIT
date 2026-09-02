/*
 * Native DSL2 stage processes.  Each process is bound to one METAHICT module
 * and declares the material files it consumes and produces.  The module
 * implementation remains the scientific reference code; Nextflow owns the
 * dependency graph, task resources, environments, and published outputs.
 */

def moduleSettings(String name) {
    def configured = params.modules
    if (configured == null) {
        return [:]
    }
    if (!(configured instanceof Map)) {
        error "The configuration key 'modules' must be a map, not ${configured.getClass().simpleName}"
    }
    def settings = configured[name]
    if (settings == null) {
        return [:]
    }
    if (!(settings instanceof Map)) {
        error "The configuration key 'modules.${name}' must be a map"
    }
    return settings
}

def modulePathSettings(String name, Collection groups = []) {
    def settings = moduleSettings(name)
    def traversed = [name]
    groups.each { group ->
        traversed << group.toString()
        def nested = settings[group]
        if (nested == null) {
            settings = [:]
        } else if (!(nested instanceof Map)) {
            error "The configuration key 'modules.${traversed.join('.')}' must be a map"
        } else {
            settings = nested
        }
    }
    return settings
}

def moduleResource(String name, String key, fallback) {
    def resources = params.resources
    if (resources instanceof Map && resources[name] instanceof Map) {
        def configured = resources[name][key]
        if (configured != null && configured.toString().trim()) {
            return configured
        }
    }
    return fallback
}

def moduleThreads(String name, Integer fallback) {
    if (params.threads != null) {
        def override = params.threads.toString().toInteger()
        if (override < 1) {
            error "--threads must be at least 1"
        }
        return override
    }
    def configured = moduleResource(name, 'threads', fallback).toString().toInteger()
    if (configured < 1) {
        error "resources.${name}.threads must be at least 1"
    }
    return configured
}

def moduleMemory(String name, String fallback) {
    if (params.memory != null) {
        def override = params.memory as nextflow.util.MemoryUnit
        if (override.toBytes() < 1) {
            error "--memory must be greater than zero"
        }
        return override
    }
    def configured = moduleResource(name, 'memory', fallback) as nextflow.util.MemoryUnit
    if (configured.toBytes() < 1) {
        error "resources.${name}.memory must be greater than zero"
    }
    return configured
}

def settingsBoolean(Map settings, String path, String key, Boolean fallback) {
    if (!settings.containsKey(key)) {
        return fallback
    }
    def value = settings[key]
    if (!(value instanceof Boolean)) {
        error "The configuration key '${path}.${key}' must be true or false"
    }
    return value
}

def modulePathBoolean(String name, Collection groups, String key, Boolean fallback) {
    def suffix = groups ? ".${groups.join('.')}" : ''
    return settingsBoolean(modulePathSettings(name, groups), "modules.${name}${suffix}", key, fallback)
}

def shellQuote(value) {
    return "'${value.toString().replace("'", "'\"'\"'")}'"
}

def settingsArguments(Map settings, Collection excluded = [], Collection falseFlags = [], Map optionAliases = [:]) {
    def ignored = excluded as Set
    def negatable = falseFlags as Set
    settings.collect { key, value ->
        def optionName = optionAliases.containsKey(key.toString()) ? optionAliases[key.toString()] : key
        def option = optionName.toString().replace('_', '-')
        if (ignored.contains(key.toString()) || value == null || value.toString() == '') {
            return ''
        }
        if (value instanceof Boolean) {
            if (value) {
                return "--${option}"
            }
            return negatable.contains(key.toString()) ? "--no-${option}" : ''
        }
        return "--${option}=${shellQuote(value)}"
    }.findAll { it }.join(' ')
}

def moduleArguments(String name, Collection excluded = [], Collection falseFlags = [], Map optionAliases = [:]) {
    return settingsArguments(moduleSettings(name), excluded, falseFlags, optionAliases)
}

def modulePathArguments(String name, Collection groups, Collection excluded = [], Collection falseFlags = [], Map optionAliases = [:]) {
    return settingsArguments(modulePathSettings(name, groups), excluded, falseFlags, optionAliases)
}

def alignmentConfigurationArguments(String name, Collection prefix = []) {
    def parts = []
    parts << modulePathArguments(name, prefix + ['bwa'], [], [], [options: 'bwa_options'])
    parts << modulePathArguments(
        name,
        prefix + ['filtering'],
        [],
        [],
        [
            min_mapping_quality: 'mapq',
            min_intra_contig_distance: 'min_intra_dist',
            min_aligned_length: 'min_match_len',
        ]
    )
    parts << modulePathArguments(
        name,
        prefix + ['sorting'],
        [],
        [],
        [memory_per_thread: 'sort_memory', temporary_directory: 'tmp_dir']
    )
    def metricsEnabled = modulePathBoolean(name, prefix + ['metrics'], 'enabled', true)
    if (!metricsEnabled) {
        parts << '--skip-metrics'
    }
    return parts.findAll { it }.join(' ')
}

def contactConfiguration(String name, Collection prefix = []) {
    def base = modulePathSettings(name, prefix)
    def path = "modules.${name}${prefix ? '.' + prefix.join('.') : ''}"
    def method = (base.method ?: 'normcc').toString().toLowerCase()
    if (!(method in ['raw', 'normcc', 'hiczin', 'bin3c', 'metator'])) {
        error "Invalid ${path}.method '${method}'"
    }
    def parts = []
    parts << modulePathArguments(
        name,
        prefix + ['raw_contacts'],
        [],
        [],
        [
            min_contact_signal: 'metacc_min_signal',
            min_contig_length: 'metacc_min_len',
            min_mapping_quality: 'metacc_min_mapq',
            min_aligned_length: 'metacc_min_match',
        ]
    )
    parts << modulePathArguments(name, prefix + ['denoising'])
    if (method != 'raw') {
        parts << modulePathArguments(
            name,
            prefix + [method],
            [],
            [],
            [max_iterations: 'max_iter', convergence_tolerance: 'tol']
        )
    }
    return [method: method, arguments: parts.findAll { it }.join(' ')]
}

def mgeAnalysisArguments() {
    def common = moduleArguments(
        'mge', ['keep_intermediates', 'genomad', 'ccfind', 'pairs', 'alignment', 'contact'], [],
        [temporary_directory: 'tmp_dir']
    )
    def genomad = modulePathArguments(
        'mge', ['genomad'], [], ['cleanup'],
        [
            splits: 'genomad_splits', sensitivity: 'genomad_sensitivity', cleanup: 'genomad_cleanup',
            restart: 'genomad_restart', preset: 'genomad_preset', min_score: 'genomad_min_score',
            max_false_discovery_rate: 'genomad_max_fdr', extra_args: 'genomad_extra_args',
        ]
    )
    def ccfind = modulePathArguments(
        'mge', ['ccfind'], [], [],
        [
            terminal_fragment_size: 'ccfind_terminal_fragment_size',
            min_percent_identity: 'ccfind_min_identity', min_aligned_length: 'ccfind_min_aligned_length',
        ]
    )
    def pairs = modulePathArguments(
        'mge', ['pairs'], [], [], [filter_method: 'pair_filter']
    )
    return [common, genomad, ccfind, pairs].findAll { it }.join(' ')
}

def mgeIntermediatePublishPath(filename, String stage) {
    if (!modulePathBoolean('mge', [], 'keep_intermediates', false)) {
        return null
    }
    def internal = "mge_${stage}"
    return filename == internal
        ? "10_MGE/intermediates/${stage}"
        : filename.replaceFirst("^${internal}/", "10_MGE/intermediates/${stage}/")
}

def restrictionEnzyme(row, String sample, String stage) {
    def sampleValue = row.enzyme?.toString()?.trim()
    if (!sampleValue) {
        error "Missing required samplesheet enzyme for ${sample}:${stage}"
    }
    return sampleValue
}

def sampleLongReadType(row) {
    def value = row.long_read_type?.toString()?.trim()
    if (!value) {
        return null
    }
    def allowed = [
        'pacbio-raw', 'pacbio-corr', 'pacbio-hifi',
        'nano-raw', 'nano-corr', 'nano-hq',
    ]
    if (!(value in allowed)) {
        error "Invalid samplesheet long_read_type '${value}'"
    }
    if (row.sg2?.toString()?.trim()) {
        error 'Samplesheet sg2 must be empty when long_read_type is set'
    }
    return value
}

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
    "${project_path}/metahict" verify
    touch conda_bundle.ok
    """
}

process PREPROCESSING {
    tag "${sample}:${library}"
    cpus { moduleThreads('preprocessing', 8) }
    memory { moduleMemory('preprocessing', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}/1_preprocessing" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'preprocessing' ? library : "${library}/${filename.replaceFirst('^preprocessing/', '')}" }

    input:
    tuple val(sample), val(row), val(library), path(read1), path(read2)

    output:
    tuple val(sample), val(row), val(library), path('preprocessing'), emit: results

    script:
    def common = moduleArguments(
        'preprocessing', ['libraries', 'trimming', 'quality_control'], [],
        [output_prefix: 'prefix']
    )
    def trimming = modulePathArguments(
        'preprocessing', ['trimming'], [], [],
        [
            min_read_length: 'minlen', quality_trim_threshold: 'trimq',
            quality_trim_ends: 'qtrim', force_trim_left: 'ftl',
            length_modulo: 'ftm', adapter_trim_end: 'ktrim',
            adapter_kmer_length: 'k', min_adapter_kmer_length: 'mink',
            adapter_hamming_distance: 'hdist', adapter_reference: 'adapter_ref',
        ]
    )
    def runPreQc = modulePathBoolean('preprocessing', ['quality_control'], 'run_before_trimming', true)
    def runPostQc = modulePathBoolean('preprocessing', ['quality_control'], 'run_after_trimming', true)
    def qc = "${runPreQc ? '' : '--skip-pre-qc-report'} ${runPostQc ? '' : '--skip-post-qc-report'}"
    def libraryKey = library == 'hic' ? 'hic' : 'shotgun'
    def useDedup = modulePathBoolean('preprocessing', ['libraries', libraryKey], 'deduplicate', library == 'hic')
    def dedup = useDedup ? '--dedup' : '--no-dedup'
    def memoryGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8).intValue())
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    mkdir -p preprocessing
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=${memoryGb} GB BBTools Java heap"
    python "${params.project_path}/modules/1_preprocessing/preprocessing.py" \\
      -p "${params.project_path}" -1 "${read1}" -2 "${read2}" \\
      -o preprocessing -t "${task.cpus}" ${dedup} ${common} ${trimming} ${qc} --xmx "${memoryGb}g"
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
    cpus { moduleThreads('assembly', 16) }
    memory { moduleMemory('assembly', '64 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'assembly' ? '2_assembly' : filename.replaceFirst('^assembly/', '2_assembly/') }

    input:
    tuple val(sample), val(row), path(shotgun_input)

    output:
    tuple val(sample), val(row), path('assembly'), emit: results

    script:
    def settings = moduleSettings('assembly')
    def longReadType = sampleLongReadType(row)
    def configuredAssembler = (settings.assembler ?: 'megahit').toString().toLowerCase()
    def assembler = longReadType ? 'metaflye' : configuredAssembler
    def assemblerOptions = [megahit: '--megahit', metaspades: '--metaspades', metaflye: '--metaflye']
    if (!assemblerOptions.containsKey(assembler)) {
        error "Invalid modules.assembly.assembler '${assembler}'; use megahit, metaspades, or metaflye"
    }
    if (!longReadType && assembler == 'metaflye') {
        error 'modules.assembly.assembler=metaflye requires a samplesheet long_read_type'
    }
    def assemblerArg = assemblerOptions[assembler]
    def common = moduleArguments(
        'assembly', ['assembler', 'quality_control', 'megahit', 'metaspades', 'metaflye'], [],
        [min_contig_length: 'min_len', temporary_directory: 'tmp_dir', keep_temporary_files: 'keep_temp']
    )
    def methodAliases = [
        min_kmer_length: 'k_min', max_kmer_length: 'k_max', kmer_step: 'k_step',
        extra_args: "${assembler}_extra_args", kmer_lengths: 'k_list',
    ]
    def methodArgs = modulePathArguments('assembly', [assembler], [], [], methodAliases)
    def runQuast = modulePathBoolean('assembly', ['quality_control'], 'run_quast', true)
    def qualityControl = runQuast ? '' : '--skip-quast'
    def memoryGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8).intValue())
    def inputSetup = longReadType
        ? "test -s \"${shotgun_input}\""
        : "sg1=\$(find -L \"${shotgun_input}\" -name 'final_*_1.fastq.gz' -print -quit)\n    sg2=\$(find -L \"${shotgun_input}\" -name 'final_*_2.fastq.gz' -print -quit)\n    test -n \"\${sg1}\" && test -n \"\${sg2}\""
    def inputArguments = longReadType
        ? "--long-reads \"${shotgun_input}\" --long-read-type ${shellQuote(longReadType)}"
        : '-1 "${sg1}" -2 "${sg2}"'
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    ${inputSetup}
    mkdir -p assembly
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=${memoryGb} GB assembler limit"
    python "${params.project_path}/modules/2_assembly/assembly.py" \\
      -p "${params.project_path}" ${inputArguments} \\
      -o assembly -t "${task.cpus}" ${assemblerArg} ${common} ${methodArgs} ${qualityControl} --memory "${memoryGb}"
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
    cpus { moduleThreads('alignment', 16) }
    memory { moduleMemory('alignment', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'alignment' ? '3_alignment' : filename.replaceFirst('^alignment/', '3_alignment/') }

    input:
    tuple val(sample), val(row), path(assembly_dir), path(hic_dir)

    output:
    tuple val(sample), val(row), path('alignment'), emit: results

    script:
    def extra = alignmentConfigurationArguments('alignment')
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -n "\${hic1}" && test -n "\${hic2}"
    mkdir -p alignment
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/3_alignment/alignment.py" \\
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
    cpus { moduleThreads('coverage', 16) }
    memory { moduleMemory('coverage', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'coverage' ? '4_coverage' : filename.replaceFirst('^coverage/', '4_coverage/') }

    input:
    tuple val(sample), val(row), path(assembly_dir), path(shotgun_input)

    output:
    tuple val(sample), val(row), path('coverage'), emit: results

    script:
    def common = moduleArguments(
        'coverage', ['mapping', 'depth'], [],
        [temporary_directory: 'tmp_dir', keep_temporary_files: 'keep_temp']
    )
    def longReadType = sampleLongReadType(row)
    def mapping = modulePathArguments(
        'coverage', ['mapping'],
        longReadType ? ['min_percent_identity'] : ['long_read_min_percent_identity'],
        [],
        [min_percent_identity: 'percent_identity', long_read_min_percent_identity: 'percent_identity']
    )
    def depth = modulePathArguments(
        'coverage', ['depth'], [], [],
        [
            min_mapping_quality: 'min_mapq', mapping_quality_weight: 'weight_mapq',
            max_excluded_edge_bases: 'max_edge_bases',
        ]
    )
    def memoryGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8).intValue())
    def inputSetup = longReadType
        ? "test -s \"${shotgun_input}\""
        : "sg1=\$(find -L \"${shotgun_input}\" -name 'final_*_1.fastq.gz' -print -quit)\n    sg2=\$(find -L \"${shotgun_input}\" -name 'final_*_2.fastq.gz' -print -quit)\n    test -n \"\${sg1}\" && test -n \"\${sg2}\""
    def inputArguments = longReadType
        ? "--long-reads \"${shotgun_input}\""
        : '-1 "${sg1}" -2 "${sg2}"'
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/metabat2/bin:\$PATH"
    ${inputSetup}
    mkdir -p coverage
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=${memoryGb} GB BBMap Java heap"
    python "${params.project_path}/modules/4_coverage/coverage.py" \\
      -p "${params.project_path}" -r "${assembly_dir}/final_assembly.fasta" \\
      ${inputArguments} -o coverage -t "${task.cpus}" ${common} ${mapping} ${depth} --memory "${memoryGb}g"
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
    cpus { moduleThreads('contact', 1) }
    memory { moduleMemory('contact', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'contact' ? '5_contact' : filename.replaceFirst('^contact/', '5_contact/') }

    input:
    tuple val(sample), val(row), path(assembly_dir), path(alignment_dir)

    output:
    tuple val(sample), val(row), path('contact'), emit: results

    script:
    def contactConfig = contactConfiguration('contact')
    def method = contactConfig.method
    def extra = contactConfig.arguments
    def enzyme = shellQuote(restrictionEnzyme(row, sample, 'contact'))
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    mkdir -p contact
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/5_contact/contact.py" ${method} \\
      -p "${params.project_path}" --bam "${alignment_dir}/sorted_map.bam" \\
      --fasta "${assembly_dir}/final_assembly.fasta" --out contact --enzyme ${enzyme} ${extra}
    test -s contact/Raw_contact_matrix.npz
    test -s contact/denoised_contact_matrix_${method}.npz
    """

    stub:
    """
    mkdir -p contact
    touch contact/Raw_contact_matrix.npz \
      contact/denoised_contact_matrix_raw.npz \
      contact/denoised_contact_matrix_normcc.npz \
      contact/denoised_contact_matrix_hiczin.npz \
      contact/denoised_contact_matrix_bin3c.npz \
      contact/denoised_contact_matrix_metator.npz \
      contact/contig_info.csv
    """
}

process BINNING {
    tag "${sample}:binning"
    cpus { moduleThreads('binning', 16) }
    memory { moduleMemory('binning', '64 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'binning' ? '6_binning' : filename.replaceFirst('^binning/', '6_binning/') }

    input:
    tuple val(sample), val(row), path(assembly_dir), path(alignment_dir)

    output:
    tuple val(sample), val(row), path('binning'), emit: results

    script:
    def common = moduleArguments(
        'binning', ['metacc', 'bin3c', 'imputecc', 'refinement', 'output'], [],
        [random_seed: 'seed', temporary_directory: 'tmp_dir', keep_temporary_files: 'keep_temp']
    )
    def metacc = modulePathArguments(
        'binning', ['metacc'], [], [],
        [
            min_contig_length: 'metacc_min_len', min_contact_signal: 'metacc_min_signal',
            min_mapping_quality: 'metacc_min_mapq', min_aligned_length: 'metacc_min_match',
            min_bin_size: 'metacc_min_binsize', normcc_discard_fraction: 'normcc_thres',
            marker_gene_count: 'num_gene',
        ]
    )
    def bin3c = modulePathArguments(
        'binning', ['bin3c'], [], [],
        [
            min_contig_length: 'bin3c_min_len', min_contact_signal: 'bin3c_min_signal',
            min_mapping_quality: 'bin3c_min_mapq', min_aligned_length: 'bin3c_min_match',
            min_cluster_extent: 'bin3c_min_extent',
        ]
    )
    def imputecc = modulePathArguments(
        'binning', ['imputecc'], [], [],
        [
            random_walk_restart_probability: 'imputecc_rwr_restart_probability',
            random_walk_threshold: 'imputecc_rwr_threshold', max_markers: 'imputecc_max_markers',
            intra_bin_threshold: 'imputecc_intra_bin_threshold', inter_bin_threshold: 'imputecc_inter_bin_threshold',
            min_bin_size: 'imputecc_min_bin_size', contamination_weight: 'imputecc_contamination_weight',
            min_completeness: 'imputecc_min_completeness', max_contamination: 'imputecc_max_contamination',
            report_quality_threshold: 'imputecc_report_quality_threshold', gene_coverage: 'imputecc_gene_coverage',
        ]
    )
    def refinementSettings = modulePathSettings('binning', ['refinement'])
    def refinement = settingsArguments(
        refinementSettings,
        ['run_checkm2', 'run_refinement', 'run_consolidation', 'ambiguous_contigs']
    )
    def refinementFlags = []
    if (!settingsBoolean(refinementSettings, 'modules.binning.refinement', 'run_checkm2', true)) refinementFlags << '--skip-checkm2'
    if (!settingsBoolean(refinementSettings, 'modules.binning.refinement', 'run_refinement', true)) refinementFlags << '--skip-refinement'
    if (!settingsBoolean(refinementSettings, 'modules.binning.refinement', 'run_consolidation', true)) refinementFlags << '--skip-consolidation'
    def ambiguous = (refinementSettings.ambiguous_contigs ?: 'best').toString().toLowerCase()
    if (!(ambiguous in ['best', 'keep', 'remove'])) {
        error "Invalid modules.binning.refinement.ambiguous_contigs '${ambiguous}'; use best, keep, or remove"
    }
    if (ambiguous == 'keep') refinementFlags << '--keep-ambiguous'
    if (ambiguous == 'remove') refinementFlags << '--remove-ambiguous'
    def outputSettings = modulePathSettings('binning', ['output'])
    def outputFlags = []
    if (!settingsBoolean(outputSettings, 'modules.binning.output', 'generate_bin3c_fasta', true)) outputFlags << '--no-fasta'
    if (!settingsBoolean(outputSettings, 'modules.binning.output', 'generate_bin3c_report', true)) outputFlags << '--no-report'
    if (!settingsBoolean(outputSettings, 'modules.binning.output', 'input_assembly_is_spades', false)) outputFlags << '--no-spades'
    if (settingsBoolean(outputSettings, 'modules.binning.output', 'only_large_bin3c_clusters', false)) outputFlags << '--only-large'
    def outputValues = settingsArguments(
        outputSettings,
        ['generate_bin3c_fasta', 'generate_bin3c_report', 'input_assembly_is_spades', 'only_large_bin3c_clusters'],
        [], [heatmap_max_image_size: 'heatmap_max_image']
    )
    def extra = ([common, metacc, bin3c, imputecc, refinement, outputValues] + refinementFlags + outputFlags).findAll { it }.join(' ')
    def enzyme = shellQuote(restrictionEnzyme(row, sample, 'binning'))
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/checkm2/bin:\$PATH"
    export CHECKM2DB="${params.checkm2_db}"
    export METAHICT_CHECKM2_BIN="${params.conda_envs_path}/checkm2/bin/checkm2"
    test -d "${params.checkm_db}"
    test -s "${params.checkm2_db}"
    mkdir -p binning
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/6_binning/run_binning.py" \\
      --modules-path "${params.project_path}/modules" \\
      --fasta "${assembly_dir}/final_assembly.fasta" --bam "${alignment_dir}/sorted_map.bam" \\
      --outdir binning --enzyme ${enzyme} -t "${task.cpus}" \\
      --checkm_db "${params.checkm_db}" ${extra}
    test -d binning/metahict/final_bins
    """

    stub:
    """
    mkdir -p binning/metahict/final_bins
    touch binning/metahict/final_bins/bin.1.fa
    printf 'bin_id\tcompleteness\tcontamination\n' > binning/metahict/final_bins_quality.tsv
    printf 'contig\tbin_id\n' > binning/metahict/contig_to_bin.tsv
    touch binning/metahict/combined_final_bins.fa
    printf 'metahict_version: 1.1.0\nmodule: binning\n' > binning/metahict/run_parameters.yaml
    """
}

process REASSEMBLY {
    tag "${sample}:reassembly"
    cpus { moduleThreads('reassembly', 16) }
    memory { moduleMemory('reassembly', '64 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'reassembly' ? '7_reassembly' : filename.replaceFirst('^reassembly/', '7_reassembly/') }

    input:
    tuple val(sample), val(row), path(assembly_dir), path(alignment_dir), path(binning_dir), path(sg_dir, stageAs: 'sg_preprocessing'), path(hic_dir, stageAs: 'hic_preprocessing')

    output:
    tuple val(sample), val(row), path('reassembly'), emit: results

    script:
    def common = moduleArguments(
        'reassembly', ['read_selection', 'read_recruitment', 'assembly', 'quality_control'], [],
        [temporary_directory: 'tmp_dir', keep_temporary_files: 'keep_temp']
    )
    def selection = modulePathArguments(
        'reassembly', ['read_selection'], [], [],
        [
            em_cutoff_quantile: 'cutoff_quantile', em_top_contigs: 'top_k',
            em_initial_n_fraction: 'em_initial_n_fraction',
            em_convergence_tolerance: 'em_convergence_tolerance',
            em_max_iterations: 'em_max_iterations',
            min_mapping_quality: 'min_mapq', min_aligned_length: 'min_match_len',
            exclude_duplicate_alignments: 'exclude_duplicates',
            write_nonselected_hic_reads: 'write_nonselected_hic',
        ]
    )
    def recruitment = modulePathArguments(
        'reassembly', ['read_recruitment'], [], [],
        [strict_mismatch_cutoff: 'strict_cut_off', permissive_mismatch_cutoff: 'permissive_cut_off']
    )
    def assemblySettings = modulePathSettings('reassembly', ['assembly'])
    def assemblyArgs = settingsArguments(
        assemblySettings, ['assemble_residual_reads'], [],
        [min_contig_length: 'min_contig_len', phred_offset: 'spades_phred_offset', extra_args: 'spades_extra_args']
    )
    if (!settingsBoolean(assemblySettings, 'modules.reassembly.assembly', 'assemble_residual_reads', true)) {
        assemblyArgs = "${assemblyArgs} --skip-residual-assembly"
    }
    def qualitySettings = modulePathSettings('reassembly', ['quality_control'])
    def quality = settingsArguments(qualitySettings, ['run_checkm2'])
    if (!settingsBoolean(qualitySettings, 'modules.reassembly.quality_control', 'run_checkm2', true)) {
        quality = "${quality} --skip-checkm2"
    }
    def extra = [common, selection, recruitment, assemblyArgs, quality].findAll { it }.join(' ')
    def memoryGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8).intValue())
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/checkm2/bin:\$PATH"
    export CHECKM2DB="${params.checkm2_db}"
    export METAHICT_CHECKM2_BIN="${params.conda_envs_path}/checkm2/bin/checkm2"
    test -s "${params.checkm2_db}"
    sg1=\$(find -L "${sg_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    sg2=\$(find -L "${sg_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    bins="${binning_dir}/metahict/final_bins"
    test -n "\${sg1}" && test -n "\${sg2}" && test -n "\${hic1}" && test -n "\${hic2}" && test -d "\${bins}"
    mkdir -p reassembly
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=${memoryGb} GB total SPAdes budget"
    python "${params.project_path}/modules/7_reassembly/run_reassembly.py" \\
      -p "${params.project_path}/modules" --bin "\${bins}" \\
      --assembly "${assembly_dir}/final_assembly.fasta" --hic1 "\${hic1}" --hic2 "\${hic2}" \\
      --sg1 "\${sg1}" --sg2 "\${sg2}" --bam "${alignment_dir}/sorted_map.bam" \\
      --outdir reassembly -t "${task.cpus}" \\
      --checkm2_db "${params.checkm2_db}" ${extra} -m "${memoryGb}"
    test -d reassembly/reassembled_bins
    """

    stub:
    """
    mkdir -p reassembly/reassembled_bins reassembly/quality reassembly/figures
    touch reassembly/reassembled_bins/bin1.fa reassembly/combined_contigs.fa
    touch reassembly/residual_contigs.fa
    printf 'bin\tcompleteness\tcontamination\n' > reassembly/reassembled_bins_quality.tsv
    printf 'final_bin\tselected_candidate\tselection_type\n' > reassembly/reassembled_bin_name_map.tsv
    printf '{"selection": {}}\n' > reassembly/read_selection_summary.json
    printf 'bin\tcompleteness\tcontamination\n' > reassembly/quality/original_bins_quality.tsv
    printf 'Name\tCompleteness\tContamination\n' > reassembly/quality/checkm2_quality_report.tsv
    printf 'metahict_version: 1.1.0\nmodule: reassembly\n' > reassembly/run_parameters.yaml
    """
}

process ANNOTATION {
    tag "${sample}:annotation"
    cpus { moduleThreads('annotation', 8) }
    memory { moduleMemory('annotation', '64 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'annotation' ? '9_annotation' : filename.replaceFirst('^annotation/', '9_annotation/') }

    input:
    tuple val(sample), val(row), path(mag_dir)

    output:
    tuple val(sample), val(row), path('annotation'), emit: results

    script:
    def settings = moduleSettings('annotation')
    def common = moduleArguments(
        'annotation', ['ani_screen', 'classification'], [],
        [
            pplacer_threads: 'pplacer_cpus', genome_extension: 'extension', output_prefix: 'prefix',
            scratch_directory: 'scratch_dir', temporary_directory: 'tmp_dir',
            continue_on_genome_error: 'force', extra_args: 'gtdbtk_extra_args',
        ]
    )
    if (settings.pplacer_threads != null && settings.pplacer_threads.toString().toInteger() > task.cpus) {
        error "modules.annotation.pplacer_threads cannot exceed the resolved annotation threads (${task.cpus})"
    }
    def aniSettings = modulePathSettings('annotation', ['ani_screen'])
    def aniEnabled = settingsBoolean(aniSettings, 'modules.annotation.ani_screen', 'enabled', false)
    def ani = aniEnabled ? '--no-skip-ani-screen' : '--skip-ani-screen'
    ani = "${ani} ${settingsArguments(aniSettings, ['enabled'], [], [mash_database: 'mash_db'])}"
    def classification = modulePathArguments(
        'annotation', ['classification'], [], [],
        [
            min_amino_acid_percent: 'min_perc_aa', min_alignment_fraction: 'min_af',
            use_full_tree: 'full_tree',
        ]
    )
    def extra = [common, ani, classification].findAll { it }.join(' ')
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/gtdbtk-2.4.0/bin:\$PATH"
    test -d "${params.gtdbtk_db}"
    mkdir -p annotation
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/9_annotation/annotation.py" \\
      -p "${params.project_path}" --mag-dir "${mag_dir}" \\
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
    tag "${sample}:${bin_id}:scaffolding"
    cpus { moduleThreads('scaffolding', 8) }
    memory { moduleMemory('scaffolding', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'scaffolding' ? "8_scaffolding/${bin_id}" : filename.replaceFirst('^scaffolding/', "8_scaffolding/${bin_id}/") }

    input:
    tuple val(sample), val(row), val(bin_id), path(bin_fasta), val(use_bam), path(bam_input), path(hic_dir)

    output:
    tuple val(sample), val(row), val(bin_id), path('scaffolding'), emit: results

    script:
    def runCheckm2 = modulePathBoolean('scaffolding', ['quality_control'], 'run_checkm2', true)
    def skipCheckm2 = params.scaffolding_skip_checkm2 || !runCheckm2
    def common = moduleArguments(
        'scaffolding', ['input_filter', 'alignment', 'contacts', 'yahs', 'heatmap', 'quality_control'], [],
        [temporary_directory: 'tmp_dir', keep_temporary_files: 'keep_temp']
    )
    def inputFilter = modulePathArguments(
        'scaffolding', ['input_filter'], [], [], [min_contig_length: 'min_contig_len']
    )
    def alignment = modulePathArguments(
        'scaffolding', ['alignment'], [], [], [sort_memory_per_thread: 'sort_memory']
    )
    def contacts = modulePathArguments(
        'scaffolding', ['contacts'], [], [],
        [
            min_mapping_quality: 'metacc_min_mapq', min_contig_length: 'metacc_min_len',
            min_aligned_length: 'metacc_min_match', min_contact_signal: 'metacc_min_signal',
            normcc_discard_fraction: 'normcc_thres',
        ]
    )
    def yahsSettings = modulePathSettings('scaffolding', ['yahs'])
    def yahs = settingsArguments(
        yahsSettings, ['contig_error_correction', 'scaffold_error_correction', 'memory_check'], [],
        [
            resolutions: 'yahs_resolutions', min_mapping_quality: 'yahs_min_mapq',
            min_contig_length: 'yahs_min_contig_len', rounds: 'yahs_rounds', extra_args: 'yahs_extra_args',
        ]
    )
    if (!settingsBoolean(yahsSettings, 'modules.scaffolding.yahs', 'contig_error_correction', true)) yahs = "${yahs} --yahs-no-contig-ec"
    if (!settingsBoolean(yahsSettings, 'modules.scaffolding.yahs', 'scaffold_error_correction', true)) yahs = "${yahs} --yahs-no-scaffold-ec"
    if (!settingsBoolean(yahsSettings, 'modules.scaffolding.yahs', 'memory_check', true)) yahs = "${yahs} --yahs-no-mem-check"
    def heatmap = modulePathArguments(
        'scaffolding', ['heatmap'], [], [],
        [segment_resolution: 'resolution', max_image_size: 'heatmap_max_image']
    )
    def extra = [common, inputFilter, alignment, contacts, yahs, heatmap].findAll { it }.join(' ')
    def enzyme = shellQuote(restrictionEnzyme(row, sample, 'scaffolding'))
    def checkm2Setup = skipCheckm2 ? '' : """
    export CHECKM2DB="${params.checkm2_db}"
    export METAHICT_CHECKM2_BIN="${params.conda_envs_path}/checkm2/bin/checkm2"
    test -s "${params.checkm2_db}"
    """
    def checkm2Args = skipCheckm2 ? '--skip-checkm2' : "--checkm2_db \"${params.checkm2_db}\""
    def bamArg = use_bam ? "--bam \"${bam_input}\"" : ''
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/checkm2/bin:\$PATH"
    ${checkm2Setup}
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -n "\${hic1}" && test -n "\${hic2}" && test -s "${bin_fasta}"
    mkdir -p scaffolding
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/8_scaffolding/scaffolding.py" \\
      -p "${params.project_path}/modules" --fasta "${bin_fasta}" ${bamArg} --enzyme ${enzyme} \\
      --outdir scaffolding --hic1 "\${hic1}" --hic2 "\${hic2}" \\
      -t "${task.cpus}" \\
      ${checkm2Args} ${extra}
    test -d scaffolding
    """

    stub:
    """
    mkdir -p scaffolding
    touch scaffolding/scaffolded_bin.fa scaffolding/scaffolded_bin.agp
    printf 'Scaffolding Quality Metrics\n' > scaffolding/scaffolding_metrics.txt
    printf 'bin\tstatus\treason\ttotal_contigs\teligible_contigs\tlongest_contig_bp\tmin_contig_length\n%s\tcompleted\t\t2\t2\t5000\t5000\n' "${bin_fasta.name}" > scaffolding/scaffolding_status.tsv
    """
}

process MGE_ALIGNMENT {
    tag "${sample}:mge-alignment"
    cpus { moduleThreads('mge', 16) }
    memory { moduleMemory('mge', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> mgeIntermediatePublishPath(filename, 'alignment') }

    input:
    tuple val(sample), val(row), path(fasta), path(hic_dir)

    output:
    tuple val(sample), val(row), path('mge_alignment'), emit: results

    script:
    def extra = alignmentConfigurationArguments('mge', ['alignment'])
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    hic1=\$(find -L "${hic_dir}" -name 'final_*_1.fastq.gz' -print -quit)
    hic2=\$(find -L "${hic_dir}" -name 'final_*_2.fastq.gz' -print -quit)
    test -s "${fasta}"
    mkdir -p mge_alignment
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/3_alignment/alignment.py" \\
      -p "${params.project_path}" -r "${fasta}" \\
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
    // Contact normalization is serial; keep the shared MGE memory allocation
    // but do not reserve idle CPUs for this internal step.
    cpus 1
    memory { moduleMemory('mge', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> mgeIntermediatePublishPath(filename, 'contact') }

    input:
    tuple val(sample), val(row), path(fasta), path(mge_alignment_dir)

    output:
    tuple val(sample), val(row), path('mge_contact'), emit: results

    script:
    def contactConfig = contactConfiguration('mge', ['contact'])
    def method = contactConfig.method
    def extra = contactConfig.arguments
    def enzyme = shellQuote(restrictionEnzyme(row, sample, 'mge-contact'))
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:\$PATH"
    mkdir -p mge_contact
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/5_contact/contact.py" ${method} \\
      -p "${params.project_path}" --bam "${mge_alignment_dir}/sorted_map.bam" \\
      --fasta "${fasta}" --out mge_contact --enzyme ${enzyme} ${extra}
    test -s mge_contact/Raw_contact_matrix.npz
    test -s mge_contact/denoised_contact_matrix_${method}.npz
    """

    stub:
    """
    mkdir -p mge_contact
    touch mge_contact/Raw_contact_matrix.npz \
      mge_contact/denoised_contact_matrix_raw.npz \
      mge_contact/denoised_contact_matrix_normcc.npz \
      mge_contact/denoised_contact_matrix_hiczin.npz \
      mge_contact/denoised_contact_matrix_bin3c.npz \
      mge_contact/denoised_contact_matrix_metator.npz \
      mge_contact/contig_info.csv
    """
}

process MGE {
    tag "${sample}:mge"
    cpus { moduleThreads('mge', 16) }
    memory { moduleMemory('mge', '32 GB') }
    conda params.metahict_conda
    publishDir { "${params.out_root}/${sample}" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'mge' ? '10_MGE' : filename.replaceFirst('^mge/', '10_MGE/') }

    input:
    tuple val(sample), val(row), path(fasta), path(host_dir), path(mge_contact_dir)

    output:
    path('mge'), emit: results

    script:
    def extra = mgeAnalysisArguments()
    def contactMethod = contactConfiguration('mge', ['contact']).method
    """
    set -euo pipefail
    export CONDA_ENVS_PATH="${params.conda_envs_path}:\${CONDA_ENVS_PATH:-}"
    export PATH="${params.conda_envs_path}/metahict_env/bin:${params.conda_envs_path}/genomad/bin:${params.conda_envs_path}/ccfind_env/bin:\$PATH"
    test -d "${params.genomad_db}"
    mkdir -p mge
    echo "[INFO] Resources: threads=${task.cpus}; task_memory=${task.memory}; tool_memory=executor allocation"
    python "${params.project_path}/modules/10_MGE/mge.py" \\
      -p "${params.project_path}" \\
      --fasta "${fasta}" --host-dir "${host_dir}" \\
      --contact "${mge_contact_dir}/denoised_contact_matrix_${contactMethod}.npz" \\
      --raw-contact "${mge_contact_dir}/Raw_contact_matrix.npz" \\
      --outdir mge -t "${task.cpus}" --genomad_db "${params.genomad_db}" \\
      ${extra}
    test -d mge
    """

    stub:
    """
    mkdir -p mge/input_assembly mge/mge_reports
    touch mge/candidate_mge_host_pairs_zscore_filtered.tsv \\
      mge/sequence_topology.tsv mge/input_assembly/contig_to_host.tsv \\
      mge/mge_reports/MGE_summary.txt
    """
}
