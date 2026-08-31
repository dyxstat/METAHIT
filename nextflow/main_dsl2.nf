#!/usr/bin/env nextflow

/*
 * METAHICT native DSL2 entry point.
 *
 * Each process invokes a numbered Python stage directly. The workflow
 * materialises every scientific stage boundary and its file contract.
 */
nextflow.enable.dsl = 2

params.samplesheet = null
params.entry_module = 'all'
params.project_path = "${baseDir}/.."
params.out_root = "${baseDir}/../nextflow_results_dsl2"
params.preprocessing_dir = null
params.assembly_dir = null
params.alignment_dir = null
params.binning_dir = null
params.annotation_mag_dir = null
params.scaffolding_bin = null
params.scaffolding_bam = null
params.scaffolding_skip_checkm2 = false
params.mge_fasta = null
params.mge_host_dir = null
params.mge_alignment_dir = null
params.mge_contact_dir = null
params.report_dir = "${baseDir}/reports_dsl2"
params.memory = null
params.metahict_conda = "${baseDir}/../installation/env.yaml"
params.conda_envs_path = "${baseDir}/../conda_envs"
params.require_conda_bundle = false

include { CONDA_BUNDLE; PREPROCESSING as PREPROCESSING_SG; PREPROCESSING as PREPROCESSING_HIC; ASSEMBLY; ALIGNMENT; COVERAGE; CONTACT; BINNING; REASSEMBLY; ANNOTATION; SCAFFOLDING; MGE_ALIGNMENT; MGE_CONTACT; MGE } from './modules/local/metahict_modules'

def normaliseRow(row) {
    row.collectEntries { key, value ->
        [(key.toString().trim()): value == null ? '' : value.toString().trim()]
    }
}

def inputFile(value, label) {
    if (!value) {
        error "Samplesheet value is required: ${label}"
    }
    def path = value.startsWith('/') ? value : "${params.project_path}/${value}"
    return file(path)
}

def rowLongReadType(row) {
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

def databasePath(value, label, requiredType) {
    if (!value) {
        error "Missing required database parameter: --${label}"
    }
    def configured = value.startsWith('/') ? value : "${params.project_path}/${value}"
    def database = file(configured)
    if (!database.exists()) {
        error "Configured ${label} path does not exist: ${database}"
    }
    if (requiredType == 'file' && !database.isFile()) {
        error "Configured ${label} must be a file: ${database}"
    }
    if (requiredType == 'directory' && !database.isDirectory()) {
        error "Configured ${label} must be a directory: ${database}"
    }
    return database.toAbsolutePath().toString()
}

def resolveDatabaseParams(required) {
    if (required.contains('checkm')) {
        params.checkm_db = databasePath(params.checkm_db, 'checkm_db', 'directory')
    }
    if (required.contains('checkm2')) {
        params.checkm2_db = databasePath(params.checkm2_db, 'checkm2_db', 'file')
    }
    if (required.contains('gtdbtk')) {
        params.gtdbtk_db = databasePath(params.gtdbtk_db, 'gtdbtk_db', 'directory')
    }
    if (required.contains('genomad')) {
        params.genomad_db = databasePath(params.genomad_db, 'genomad_db', 'directory')
    }
}

def samplesheetRows() {
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }
    Channel
        .fromPath(params.samplesheet)
        .ifEmpty { error "Samplesheet not found: ${params.samplesheet}" }
        .splitCsv(header: true, sep: ',', quote: '"')
        .map { raw ->
            def row = normaliseRow(raw)
            if (!row.sample) {
                error 'Every samplesheet row requires a non-empty sample column'
            }
            tuple(row.sample, row)
        }
}

def preprocessingLibraryBoolean(String library, String key, Boolean fallback) {
    def modules = params.modules
    if (modules == null || !(modules instanceof Map) || modules.preprocessing == null) {
        return fallback
    }
    if (!(modules.preprocessing instanceof Map)) {
        error "The configuration key 'modules.preprocessing' must be a map"
    }
    def libraries = modules.preprocessing.libraries
    if (libraries == null) {
        return fallback
    }
    if (!(libraries instanceof Map) || !(libraries[library] instanceof Map)) {
        error "The configuration key 'modules.preprocessing.libraries.${library}' must be a map"
    }
    if (!libraries[library].containsKey(key)) {
        return fallback
    }
    def value = libraries[library][key]
    if (!(value instanceof Boolean)) {
        error "The configuration key 'modules.preprocessing.libraries.${library}.${key}' must be true or false"
    }
    return value
}

def selectedPreprocessingLibraries() {
    def selected = []
    if (preprocessingLibraryBoolean('shotgun', 'enabled', true)) {
        selected << 'sg'
    }
    if (preprocessingLibraryBoolean('hic', 'enabled', true)) {
        selected << 'hic'
    }
    if (!selected) {
        error "At least one preprocessing library must be enabled under modules.preprocessing.libraries"
    }
    return selected
}

def preprocessingLibraries() {
    def selected = selectedPreprocessingLibraries()
    samplesheetRows()
        .flatMap { sample, row ->
            selected
                .findAll { library -> !(library == 'sg' && rowLongReadType(row)) }
                .collect { library ->
                def read1_key = library == 'sg' ? 'sg1' : 'hic1'
                def read2_key = library == 'sg' ? 'sg2' : 'hic2'
                tuple(
                    sample,
                    row,
                    library,
                    inputFile(row[read1_key], "${sample}:${read1_key}"),
                    inputFile(row[read2_key], "${sample}:${read2_key}")
                )
            }
        }
}

def longReadShotgunChannel() {
    samplesheetRows()
        .filter { sample, row -> rowLongReadType(row) }
        .map { sample, row ->
            tuple(sample, row, inputFile(row.sg1, "${sample}:sg1"))
        }
}

def shortReadShotgunStageChannel() {
    samplesheetRows()
        .filter { sample, row -> !rowLongReadType(row) }
        .map { sample, row ->
            if (!params.preprocessing_dir) {
                error 'Short-read assembly and coverage require --preprocessing_dir PATH'
            }
            def configured = "${params.preprocessing_dir}/sg".replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            tuple(sample, row, file(resolved))
        }
}

def shotgunStageChannel() {
    shortReadShotgunStageChannel().mix(longReadShotgunChannel())
}

def resolveStageDirectory(String configuredPath, markerSpec, legacySubdirectories = []) {
    def markers = markerSpec instanceof Collection ? markerSpec : [markerSpec]
    def base = file(configuredPath)
    def candidates = [base]
    legacySubdirectories.each { legacy -> candidates << file("${configuredPath}/${legacy}") }
    def selected = candidates.find { candidate ->
        markers.any { marker -> file("${candidate}/${marker}").exists() }
    }
    if (!selected) {
        def checked = candidates.collectMany { candidate ->
            markers.collect { marker -> "${candidate}/${marker}" }
        }.join(', ')
        error "Stage directory does not contain any required marker. Checked: ${checked}"
    }
    return selected
}

def directStageChannel(String path, marker = null, legacySubdirectories = []) {
    samplesheetRows()
        .map { sample, row ->
            def configured = path.toString().replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            def directory = marker
                ? resolveStageDirectory(resolved, marker, legacySubdirectories)
                : file(resolved)
            tuple(sample, row, directory)
        }
}

def libraryStageChannel(String library) {
    if (!params.preprocessing_dir) {
        error "This selected module requires --preprocessing_dir PATH"
    }
    // Select the library directory itself. This supports both the current flat
    // layout and legacy outputs that retain a nested preprocessing/ folder,
    // because downstream stages discover the paired FASTQs recursively.
    return directStageChannel("${params.preprocessing_dir}/${library}")
}

def assemblyStageChannel() {
    if (!params.assembly_dir) {
        error "This selected module requires --assembly_dir PATH"
    }
    return directStageChannel(params.assembly_dir, 'final_assembly.fasta', ['assembly'])
}

def alignmentStageChannel() {
    if (!params.alignment_dir) {
        error "This selected module requires --alignment_dir PATH"
    }
    return directStageChannel(params.alignment_dir, 'sorted_map.bam', ['alignment'])
}

def binningStageChannel() {
    if (!params.binning_dir) {
        error "This selected module requires --binning_dir PATH"
    }
    return directStageChannel(params.binning_dir, 'metahict/final_bins', ['binning'])
}

def annotationMagStageChannel() {
    if (!params.annotation_mag_dir) {
        error "The annotation entry requires --annotation_mag_dir PATH"
    }
    samplesheetRows()
        .map { sample, row ->
            def configured = params.annotation_mag_dir.toString().replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            def mag_directory = file(resolved)
            if (!mag_directory.exists() || !mag_directory.isDirectory()) {
                error "Annotation MAG directory does not exist or is not a directory: ${mag_directory}"
            }
            tuple(sample, row, mag_directory)
        }
}

def scaffoldingBinId(bin) {
    bin.name.replaceFirst(/(?i)\.(fa|fasta|fna)$/, '')
}

def scaffoldingBinStageChannel() {
    if (!params.scaffolding_bin) {
        error "The scaffolding entry requires --scaffolding_bin PATH"
    }
    samplesheetRows()
        .map { sample, row ->
            def configured = params.scaffolding_bin.toString().replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            def bin = file(resolved)
            if (!bin.exists() || !bin.isFile()) {
                error "Scaffolding bin FASTA does not exist or is not a file: ${bin}"
            }
            def use_bam = params.scaffolding_bam != null
            def bam = file("${baseDir}/assets/NO_SCAFFOLDING_BAM")
            if (use_bam) {
                def configured_bam = params.scaffolding_bam.toString().replace('{sample}', sample.toString())
                def resolved_bam = configured_bam.startsWith('/') ? configured_bam : "${params.project_path}/${configured_bam}"
                bam = file(resolved_bam)
                if (!bam.exists() || !bam.isFile()) {
                    error "Scaffolding BAM does not exist or is not a file: ${bam}"
                }
            }
            tuple(sample, row, scaffoldingBinId(bin), bin, use_bam, bam)
        }
}

def reassembledBinStageChannel(results) {
    results.flatMap { sample, row, reassembly_directory ->
        def bins = files(
            "${reassembly_directory}/reassembled_bins/*.{fa,fasta,fna}",
            checkIfExists: true
        ) as List
        bins.sort { left, right -> left.name <=> right.name }
            .collect { bin ->
                tuple(
                    sample,
                    row,
                    scaffoldingBinId(bin),
                    bin,
                    false,
                    file("${baseDir}/assets/NO_SCAFFOLDING_BAM")
                )
            }
    }
}

def binningBinStageChannel(results) {
    results.flatMap { sample, row, binning_directory ->
        def bins = files(
            "${binning_directory}/metahict/final_bins/*.{fa,fasta,fna}",
            checkIfExists: true
        ) as List
        bins.sort { left, right -> left.name <=> right.name }
            .collect { bin ->
                tuple(
                    sample,
                    row,
                    scaffoldingBinId(bin),
                    bin,
                    false,
                    file("${baseDir}/assets/NO_SCAFFOLDING_BAM")
                )
            }
    }
}

def mgeAlignmentStageChannel() {
    if (params.mge_alignment_dir) {
        return directStageChannel(
            params.mge_alignment_dir,
            'sorted_map.bam',
            ['alignment', 'mge_alignment', 'mge_alignment/mge_alignment']
        )
    }
    return null
}

def mgeFastaStageChannel() {
    if (!params.mge_fasta) {
        error "The MGE entry requires --mge_fasta PATH"
    }
    samplesheetRows()
        .map { sample, row ->
            def configured = params.mge_fasta.toString().replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            def fasta = file(resolved)
            if (!fasta.exists() || !fasta.isFile()) {
                error "MGE FASTA does not exist or is not a file: ${fasta}"
            }
            tuple(sample, row, fasta)
        }
}

def mgeHostStageChannel() {
    if (!params.mge_host_dir) {
        error "The MGE entry requires --mge_host_dir PATH"
    }
    samplesheetRows()
        .map { sample, row ->
            def configured = params.mge_host_dir.toString().replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            def host_directory = file(resolved)
            if (!host_directory.exists() || !host_directory.isDirectory()) {
                error "MGE host directory does not exist or is not a directory: ${host_directory}"
            }
            tuple(sample, row, host_directory)
        }
}

def mgeContactStageChannel() {
    if (params.mge_contact_dir) {
        def configured = params.modules instanceof Map && params.modules.mge instanceof Map &&
            params.modules.mge.contact instanceof Map ? params.modules.mge.contact.method : null
        def method = (configured ?: 'normcc').toString().toLowerCase()
        return directStageChannel(
            params.mge_contact_dir,
            "denoised_contact_matrix_${method}.npz",
            ['contact', 'mge_contact', 'mge_contact/mge_contact']
        )
    }
    return null
}

workflow RUN_ALL {
    // Resolve database locations once, before any analytical task starts.
    // These are normal Nextflow inputs; users do not need to export database
    // environment variables.
    resolveDatabaseParams(['checkm', 'checkm2', 'gtdbtk', 'genomad'])

    def selected_libraries = selectedPreprocessingLibraries()
    if (!(selected_libraries.contains('sg') && selected_libraries.contains('hic'))) {
        error "Full workflow requires modules.preprocessing.libraries.shotgun.enabled and hic.enabled to both be true."
    }

    /* Run SG preprocessing first, then Hi-C preprocessing.  This prevents
     * large SG and Hi-C libraries from competing for memory and disk I/O. */
    libraries = preprocessingLibraries()

    if (params.require_conda_bundle) {
        CONDA_BUNDLE(Channel.value(params.project_path))
        conda_bundle_gate = CONDA_BUNDLE.out.ready
    } else {
        conda_bundle_gate = Channel.value('not-required')
    }
    gated_libraries = libraries.combine(conda_bundle_gate)
        .map { sample, row, library, read1, read2, bundle_marker -> tuple(sample, row, library, read1, read2) }
    long_shotgun = longReadShotgunChannel()
        .combine(conda_bundle_gate)
        .map { sample, row, long_reads, bundle_marker -> tuple(sample, row, long_reads) }

    sg_libraries = gated_libraries.filter { sample, row, library, read1, read2 -> library == 'sg' }
    hic_libraries = gated_libraries.filter { sample, row, library, read1, read2 -> library == 'hic' }

    PREPROCESSING_SG(sg_libraries)
    short_hic_libraries = hic_libraries
        .filter { sample, row, library, read1, read2 -> !rowLongReadType(row) }
    long_hic_libraries = hic_libraries
        .filter { sample, row, library, read1, read2 -> rowLongReadType(row) }
    short_hic_after_sg = short_hic_libraries
        .combine(PREPROCESSING_SG.out.results.collect())
        .map { values -> tuple(values[0], values[1], values[2], values[3], values[4]) }
    hic_after_sg = long_hic_libraries.mix(short_hic_after_sg)
    PREPROCESSING_HIC(hic_after_sg)

    sg_preprocessed = PREPROCESSING_SG.out.results
        .map { sample, row, library, directory -> tuple(sample, row, directory) }
    hic_preprocessed = PREPROCESSING_HIC.out.results
        .map { sample, row, library, directory -> tuple(sample, row, directory) }

    shotgun_sources = sg_preprocessed.mix(long_shotgun)
    ASSEMBLY(shotgun_sources)

    alignment_input = ASSEMBLY.out.results
        .join(hic_preprocessed)
        .map { sample, assembly_row, assembly_directory, hic_row, hic_directory ->
            tuple(sample, assembly_row, assembly_directory, hic_directory)
        }
    ALIGNMENT(alignment_input)

    coverage_input = ASSEMBLY.out.results
        .join(shotgun_sources)
        .map { sample, assembly_row, assembly_directory, shotgun_row, shotgun_input ->
            tuple(sample, assembly_row, assembly_directory, shotgun_input)
        }
    COVERAGE(coverage_input)

    contact_input = ASSEMBLY.out.results
        .join(ALIGNMENT.out.results)
        .map { sample, assembly_row, assembly_directory, alignment_row, alignment_directory ->
            tuple(sample, assembly_row, assembly_directory, alignment_directory)
        }
    CONTACT(contact_input)

    binning_input = ASSEMBLY.out.results
        .join(ALIGNMENT.out.results)
        .map { sample, assembly_row, assembly_directory, alignment_row, alignment_directory ->
            tuple(sample, assembly_row, assembly_directory, alignment_directory)
        }
    BINNING(binning_input)

    reassembly_input = BINNING.out.results
        .join(ASSEMBLY.out.results)
        .join(ALIGNMENT.out.results)
        .join(sg_preprocessed)
        .join(hic_preprocessed)
        .map { sample, binning_row, binning_directory, assembly_row, assembly_directory, alignment_row, alignment_directory, sg_row, sg_directory, hic_row, hic_directory ->
            tuple(sample, binning_row, assembly_directory, alignment_directory, binning_directory, sg_directory, hic_directory)
        }
    REASSEMBLY(reassembly_input)

    short_annotation_input = REASSEMBLY.out.results
        .map { sample, reassembly_row, reassembly_directory ->
            tuple(
                sample,
                reassembly_row,
                file("${reassembly_directory}/reassembled_bins", checkIfExists: true)
            )
        }
    long_annotation_input = BINNING.out.results
        .filter { sample, row, binning_directory -> rowLongReadType(row) }
        .map { sample, row, binning_directory ->
            tuple(
                sample,
                row,
                file("${binning_directory}/metahict/final_bins", checkIfExists: true)
            )
        }
    annotation_input = short_annotation_input.mix(long_annotation_input)
    ANNOTATION(annotation_input)

    short_scaffolding_bins = reassembledBinStageChannel(REASSEMBLY.out.results)
    long_scaffolding_bins = binningBinStageChannel(
        BINNING.out.results.filter { sample, row, binning_directory -> rowLongReadType(row) }
    )
    scaffolding_input = short_scaffolding_bins
        .mix(long_scaffolding_bins)
        .combine(hic_preprocessed, by: 0)
        .map { sample, reassembly_row, bin_id, bin_fasta, use_bam, bam_input, hic_row, hic_directory ->
            tuple(sample, reassembly_row, bin_id, bin_fasta, use_bam, bam_input, hic_directory)
        }
    SCAFFOLDING(scaffolding_input)

    /* Short-read runs use reassembled outputs. Long-read runs skip reassembly
     * and use the metaFlye assembly plus the final binning MAGs. */
    short_mge_sources = REASSEMBLY.out.results
        .map { sample, reassembly_row, reassembly_directory ->
            tuple(
                sample,
                reassembly_row,
                file("${reassembly_directory}/combined_contigs.fa", checkIfExists: true),
                file("${reassembly_directory}/reassembled_bins", checkIfExists: true)
            )
        }
    long_mge_sources = ASSEMBLY.out.results
        .join(BINNING.out.results.filter { sample, row, binning_directory -> rowLongReadType(row) })
        .map { sample, assembly_row, assembly_directory, binning_row, binning_directory ->
            tuple(
                sample,
                assembly_row,
                file("${assembly_directory}/final_assembly.fasta", checkIfExists: true),
                file("${binning_directory}/metahict/final_bins", checkIfExists: true)
            )
        }
    mge_sources = short_mge_sources.mix(long_mge_sources)
    mge_alignment_input = mge_sources
        .join(hic_preprocessed)
        .map { sample, reassembly_row, assembly_fasta, bin_directory, hic_row, hic_directory ->
            tuple(sample, reassembly_row, assembly_fasta, hic_directory)
        }
    MGE_ALIGNMENT(mge_alignment_input)

    mge_contact_input = mge_sources
        .join(MGE_ALIGNMENT.out.results)
        .map { sample, reassembly_row, assembly_fasta, bin_directory, alignment_row, alignment_directory ->
            tuple(sample, reassembly_row, assembly_fasta, alignment_directory)
        }
    MGE_CONTACT(mge_contact_input)

    mge_input = mge_sources
        .join(MGE_CONTACT.out.results)
        .map { sample, reassembly_row, assembly_fasta, bin_directory, contact_row, contact_directory ->
            tuple(sample, reassembly_row, assembly_fasta, bin_directory, contact_directory)
        }
    MGE(mge_input)
}

workflow PREPROCESSING_WORKFLOW {
    def selected_libraries = selectedPreprocessingLibraries()
    libraries = preprocessingLibraries()

    if (params.require_conda_bundle) {
        CONDA_BUNDLE(Channel.value(params.project_path))
        conda_bundle_gate = CONDA_BUNDLE.out.ready
    } else {
        conda_bundle_gate = Channel.value('not-required')
    }
    gated_libraries = libraries.combine(conda_bundle_gate)
        .map { sample, row, library, read1, read2, bundle_marker -> tuple(sample, row, library, read1, read2) }

    if (selected_libraries.contains('sg')) {
        sg_libraries = gated_libraries.filter { sample, row, library, read1, read2 -> library == 'sg' }
        PREPROCESSING_SG(sg_libraries)
    }

    if (selected_libraries.contains('hic')) {
        hic_libraries = gated_libraries.filter { sample, row, library, read1, read2 -> library == 'hic' }
        if (selected_libraries.contains('sg')) {
            short_hic_libraries = hic_libraries
                .filter { sample, row, library, read1, read2 -> !rowLongReadType(row) }
            long_hic_libraries = hic_libraries
                .filter { sample, row, library, read1, read2 -> rowLongReadType(row) }
            short_hic_after_sg = short_hic_libraries
                .combine(PREPROCESSING_SG.out.results.collect())
                .map { values -> tuple(values[0], values[1], values[2], values[3], values[4]) }
            hic_libraries = long_hic_libraries.mix(short_hic_after_sg)
        }
        PREPROCESSING_HIC(hic_libraries)
    }
}

workflow ASSEMBLY_WORKFLOW {
    ASSEMBLY(shotgunStageChannel())
}

workflow ALIGNMENT_WORKFLOW {
    alignment_input = assemblyStageChannel()
        .join(libraryStageChannel('hic'))
        .map { sample, assembly_row, assembly_directory, hic_row, hic_directory ->
            tuple(sample, assembly_row, assembly_directory, hic_directory)
        }
    ALIGNMENT(alignment_input)
}

workflow COVERAGE_WORKFLOW {
    coverage_input = assemblyStageChannel()
        .join(shotgunStageChannel())
        .map { sample, assembly_row, assembly_directory, shotgun_row, shotgun_input ->
            tuple(sample, assembly_row, assembly_directory, shotgun_input)
        }
    COVERAGE(coverage_input)
}

workflow CONTACT_WORKFLOW {
    contact_input = assemblyStageChannel()
        .join(alignmentStageChannel())
        .map { sample, assembly_row, assembly_directory, alignment_row, alignment_directory ->
            tuple(sample, assembly_row, assembly_directory, alignment_directory)
        }
    CONTACT(contact_input)
}

workflow BINNING_WORKFLOW {
    resolveDatabaseParams(['checkm', 'checkm2'])
    binning_input = assemblyStageChannel()
        .join(alignmentStageChannel())
        .map { sample, assembly_row, assembly_directory, alignment_row, alignment_directory ->
            tuple(sample, assembly_row, assembly_directory, alignment_directory)
        }
    BINNING(binning_input)
}

workflow REASSEMBLY_WORKFLOW {
    resolveDatabaseParams(['checkm2'])
    reassembly_input = binningStageChannel()
        .join(assemblyStageChannel())
        .join(alignmentStageChannel())
        .join(libraryStageChannel('sg'))
        .join(libraryStageChannel('hic'))
        .map { sample, binning_row, binning_directory, assembly_row, assembly_directory, alignment_row, alignment_directory, sg_row, sg_directory, hic_row, hic_directory ->
            tuple(sample, binning_row, assembly_directory, alignment_directory, binning_directory, sg_directory, hic_directory)
        }
    REASSEMBLY(reassembly_input)
}

workflow SCAFFOLDING_WORKFLOW {
    if (!params.scaffolding_skip_checkm2) {
        resolveDatabaseParams(['checkm2'])
    }
    scaffolding_input = scaffoldingBinStageChannel()
        .combine(libraryStageChannel('hic'), by: 0)
        .map { sample, scaffolding_row, bin_id, bin_fasta, use_bam, bam_input, hic_row, hic_directory ->
            tuple(sample, scaffolding_row, bin_id, bin_fasta, use_bam, bam_input, hic_directory)
        }
    SCAFFOLDING(scaffolding_input)
}

workflow ANNOTATION_WORKFLOW {
    resolveDatabaseParams(['gtdbtk'])
    ANNOTATION(annotationMagStageChannel())
}

workflow MGE_WORKFLOW {
    resolveDatabaseParams(['genomad'])
    mge_sources = mgeFastaStageChannel()
            .join(mgeHostStageChannel())
            .map { sample, fasta_row, fasta, host_row, host_directory ->
                tuple(sample, fasta_row, fasta, host_directory)
            }

    if (params.mge_contact_dir) {
        mge_input = mge_sources
            .join(mgeContactStageChannel())
            .map { sample, source_row, fasta, host_directory, contact_row, contact_directory ->
                tuple(sample, source_row, fasta, host_directory, contact_directory)
            }
    } else {
        if (params.mge_alignment_dir) {
            mge_alignment_results = mgeAlignmentStageChannel()
        } else {
            mge_alignment_input = mgeFastaStageChannel()
                .join(libraryStageChannel('hic'))
                .map { sample, fasta_row, fasta, hic_row, hic_directory ->
                    tuple(sample, fasta_row, fasta, hic_directory)
                }
            MGE_ALIGNMENT(mge_alignment_input)
            mge_alignment_results = MGE_ALIGNMENT.out.results
        }

        mge_contact_input = mgeFastaStageChannel()
            .join(mge_alignment_results)
            .map { sample, fasta_row, fasta, alignment_row, alignment_directory ->
                tuple(sample, fasta_row, fasta, alignment_directory)
            }
        MGE_CONTACT(mge_contact_input)

        mge_input = mge_sources
            .join(MGE_CONTACT.out.results)
            .map { sample, source_row, fasta, host_directory, contact_row, contact_directory ->
                tuple(sample, source_row, fasta, host_directory, contact_directory)
            }
    }
    MGE(mge_input)
}

workflow {
    def selected = (params.entry_module ?: 'all').toString().trim().toLowerCase()

    if (selected == 'all') {
        RUN_ALL()
    } else if (selected == 'preprocessing') {
        PREPROCESSING_WORKFLOW()
    } else if (selected == 'assembly') {
        ASSEMBLY_WORKFLOW()
    } else if (selected == 'alignment') {
        ALIGNMENT_WORKFLOW()
    } else if (selected == 'coverage') {
        COVERAGE_WORKFLOW()
    } else if (selected == 'contact') {
        CONTACT_WORKFLOW()
    } else if (selected == 'binning') {
        BINNING_WORKFLOW()
    } else if (selected == 'reassembly') {
        REASSEMBLY_WORKFLOW()
    } else if (selected == 'scaffolding') {
        SCAFFOLDING_WORKFLOW()
    } else if (selected == 'annotation') {
        ANNOTATION_WORKFLOW()
    } else if (selected == 'mge') {
        MGE_WORKFLOW()
    } else {
        error "Unknown --entry_module '${params.entry_module}'. Use all, preprocessing, assembly, alignment, coverage, contact, binning, reassembly, scaffolding, annotation, or mge."
    }
}
