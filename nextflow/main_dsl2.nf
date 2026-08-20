#!/usr/bin/env nextflow

/*
 * METAHICT native DSL2 entry point.
 *
 * The module scripts remain the implementation of each analysis.  Unlike the
 * retired wrapper (`main_wrapper_legacy.nf`), this workflow explicitly
 * materialises every stage boundary and does not use metahict_nextflow_driver.py.
 */
nextflow.enable.dsl = 2

params.samplesheet = null
params.entry_module = 'all'
params.project_path = "${baseDir}/.."
params.out_root = "${baseDir}/../nextflow_results_dsl2"
params.module_input_root = null
params.preprocessing_root = null
params.assembly_root = null
params.alignment_root = null
params.coverage_root = null
params.binning_root = null
params.reassembly_root = null
params.annotation_root = null
params.sg_preprocessing_dir = null
params.hic_preprocessing_dir = null
params.assembly_dir = null
params.alignment_dir = null
params.coverage_dir = null
params.binning_dir = null
params.reassembly_dir = null
params.mge_alignment_dir = null
params.mge_contact_dir = null
params.preprocessing_libraries = 'sg,hic'
params.report_dir = "${baseDir}/reports_dsl2"
params.metahict_conda = "${baseDir}/../installation/env.yaml"
params.conda_envs_path = "${baseDir}/../conda_envs"
params.require_conda_bundle = false
params.container_image_override = null
params.container_images = [
    all: 'metahict-all:local',
    preprocessing: 'metahict-module1:local',
    assembly: 'metahict-module2:local',
    alignment: 'metahict-module3:local',
    coverage: 'metahict-module4:local',
    contact: 'metahict-module5:local',
    binning: 'metahict-module6:local',
    reassembly: 'metahict-module7:local',
    scaffolding: 'metahict-module8:local',
    annotation: 'metahict-module9:local',
    mge: 'metahict-module10:local'
]

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

def resolveDatabaseParams() {
    params.checkm_db = databasePath(params.checkm_db, 'checkm_db', 'directory')
    params.checkm2_db = databasePath(params.checkm2_db, 'checkm2_db', 'file')
    params.gtdbtk_db = databasePath(params.gtdbtk_db, 'gtdbtk_db', 'directory')
    params.genomad_db = databasePath(params.genomad_db, 'genomad_db', 'directory')
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

def samplesWithReads() {
    samplesheetRows()
        .map { sample, row ->
            tuple(
                sample,
                row,
                inputFile(row.sg1, "${sample}:sg1"),
                inputFile(row.sg2, "${sample}:sg2"),
                inputFile(row.hic1, "${sample}:hic1"),
                inputFile(row.hic2, "${sample}:hic2")
            )
        }
}

def selectedPreprocessingLibraries() {
    def selected = (params.preprocessing_libraries ?: 'sg,hic')
        .toString()
        .split(',')
        .collect { it.trim().toLowerCase() }
        .findAll { it }
    if (!selected) {
        error "At least one preprocessing library is required: --preprocessing_libraries sg,hic|sg|hic"
    }
    selected.each { library ->
        if (!(library in ['sg', 'hic'])) {
            error "Invalid --preprocessing_libraries value '${library}'. Use sg, hic, or sg,hic."
        }
    }
    return selected.unique()
}

def preprocessingLibraries() {
    def selected = selectedPreprocessingLibraries()
    samplesheetRows()
        .flatMap { sample, row ->
            selected.collect { library ->
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

def moduleInputRoot() {
    return params.module_input_root ?: params.out_root
}

def moduleRoot(value) {
    return value ?: moduleInputRoot()
}

def stageChannel(String root, String relativePath) {
    samplesheetRows()
        .map { sample, row ->
            tuple(sample, row, file("${root}/${sample}/${relativePath}"))
        }
}

def directStageChannel(String path) {
    samplesheetRows()
        .map { sample, row ->
            def configured = path.toString().replace('{sample}', sample.toString())
            def resolved = configured.startsWith('/') ? configured : "${params.project_path}/${configured}"
            tuple(sample, row, file(resolved))
        }
}

def libraryStageChannel(String library) {
    if (library == 'sg' && params.sg_preprocessing_dir) {
        return directStageChannel(params.sg_preprocessing_dir)
    }
    if (library == 'hic' && params.hic_preprocessing_dir) {
        return directStageChannel(params.hic_preprocessing_dir)
    }
    stageChannel(moduleRoot(params.preprocessing_root), "1_preprocessing/${library}/preprocessing")
}

def assemblyStageChannel() {
    if (params.assembly_dir) {
        return directStageChannel(params.assembly_dir)
    }
    stageChannel(moduleRoot(params.assembly_root), '2_assembly/assembly')
}

def alignmentStageChannel() {
    if (params.alignment_dir) {
        return directStageChannel(params.alignment_dir)
    }
    stageChannel(moduleRoot(params.alignment_root), '3_alignment/alignment')
}

def coverageStageChannel() {
    if (params.coverage_dir) {
        return directStageChannel(params.coverage_dir)
    }
    stageChannel(moduleRoot(params.coverage_root), '4_coverage/coverage')
}

def binningStageChannel() {
    if (params.binning_dir) {
        return directStageChannel(params.binning_dir)
    }
    stageChannel(moduleRoot(params.binning_root), '6_binning/binning')
}

def reassemblyStageChannel() {
    if (params.reassembly_dir) {
        return directStageChannel(params.reassembly_dir)
    }
    stageChannel(moduleRoot(params.reassembly_root), '7_reassembly/reassembly')
}

def annotationStageChannel() {
    stageChannel(moduleRoot(params.annotation_root), '9_annotation/annotation')
}

def mgeAlignmentStageChannel() {
    if (params.mge_alignment_dir) {
        return directStageChannel(params.mge_alignment_dir)
    }
    return null
}

def mgeContactStageChannel() {
    if (params.mge_contact_dir) {
        return directStageChannel(params.mge_contact_dir)
    }
    return null
}

workflow RUN_ALL {
    // Resolve database locations once, before any analytical task starts.
    // These are normal Nextflow inputs shared by the Conda, Docker, and
    // Apptainer profiles; users never need to export container variables.
    resolveDatabaseParams()

    def selected_libraries = selectedPreprocessingLibraries()
    if (!(selected_libraries.contains('sg') && selected_libraries.contains('hic'))) {
        error "Full workflow requires both SG and Hi-C preprocessing libraries. Use --preprocessing_libraries sg,hic."
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

    sg_libraries = gated_libraries.filter { sample, row, library, read1, read2 -> library == 'sg' }
    hic_libraries = gated_libraries.filter { sample, row, library, read1, read2 -> library == 'hic' }

    PREPROCESSING_SG(sg_libraries)
    hic_after_sg = hic_libraries
        .combine(PREPROCESSING_SG.out.results.collect())
        .map { values -> tuple(values[0], values[1], values[2], values[3], values[4]) }
    PREPROCESSING_HIC(hic_after_sg)

    sg_preprocessed = PREPROCESSING_SG.out.results
        .map { sample, row, library, directory -> tuple(sample, row, directory) }
    hic_preprocessed = PREPROCESSING_HIC.out.results
        .map { sample, row, library, directory -> tuple(sample, row, directory) }

    ASSEMBLY(sg_preprocessed)

    alignment_input = ASSEMBLY.out.results
        .join(hic_preprocessed)
        .map { sample, assembly_row, assembly_directory, hic_row, hic_directory ->
            tuple(sample, assembly_row, assembly_directory, hic_directory)
        }
    ALIGNMENT(alignment_input)

    coverage_input = ASSEMBLY.out.results
        .join(sg_preprocessed)
        .map { sample, assembly_row, assembly_directory, sg_row, sg_directory ->
            tuple(sample, assembly_row, assembly_directory, sg_directory)
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

    ANNOTATION(REASSEMBLY.out.results)

    scaffolding_input = REASSEMBLY.out.results
        .join(ALIGNMENT.out.results)
        .join(hic_preprocessed)
        .map { sample, reassembly_row, reassembly_directory, alignment_row, alignment_directory, hic_row, hic_directory ->
            tuple(sample, reassembly_row, reassembly_directory, alignment_directory, hic_directory)
        }
    SCAFFOLDING(scaffolding_input)

    /* MGE analysis deliberately uses reassembled contigs, not Module 3's
     * original-assembly BAM/contact. */
    mge_alignment_input = REASSEMBLY.out.results
        .join(hic_preprocessed)
        .map { sample, reassembly_row, reassembly_directory, hic_row, hic_directory ->
            tuple(sample, reassembly_row, reassembly_directory, hic_directory)
        }
    MGE_ALIGNMENT(mge_alignment_input)

    mge_contact_input = REASSEMBLY.out.results
        .join(MGE_ALIGNMENT.out.results)
        .map { sample, reassembly_row, reassembly_directory, alignment_row, alignment_directory ->
            tuple(sample, reassembly_row, reassembly_directory, alignment_directory)
        }
    MGE_CONTACT(mge_contact_input)

    mge_input = REASSEMBLY.out.results
        .join(MGE_CONTACT.out.results)
        .map { sample, reassembly_row, reassembly_directory, contact_row, contact_directory ->
            tuple(sample, reassembly_row, reassembly_directory, contact_directory)
        }
    MGE(mge_input)
}

workflow MODULE1_PREPROCESSING {
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
            hic_libraries = hic_libraries
                .combine(PREPROCESSING_SG.out.results.collect())
                .map { values -> tuple(values[0], values[1], values[2], values[3], values[4]) }
        }
        PREPROCESSING_HIC(hic_libraries)
    }
}

workflow MODULE2_ASSEMBLY {
    ASSEMBLY(libraryStageChannel('sg'))
}

workflow MODULE3_ALIGNMENT {
    alignment_input = assemblyStageChannel()
        .join(libraryStageChannel('hic'))
        .map { sample, assembly_row, assembly_directory, hic_row, hic_directory ->
            tuple(sample, assembly_row, assembly_directory, hic_directory)
        }
    ALIGNMENT(alignment_input)
}

workflow MODULE4_COVERAGE {
    coverage_input = assemblyStageChannel()
        .join(libraryStageChannel('sg'))
        .map { sample, assembly_row, assembly_directory, sg_row, sg_directory ->
            tuple(sample, assembly_row, assembly_directory, sg_directory)
        }
    COVERAGE(coverage_input)
}

workflow MODULE5_CONTACT {
    contact_input = assemblyStageChannel()
        .join(alignmentStageChannel())
        .map { sample, assembly_row, assembly_directory, alignment_row, alignment_directory ->
            tuple(sample, assembly_row, assembly_directory, alignment_directory)
        }
    CONTACT(contact_input)
}

workflow MODULE6_BINNING {
    resolveDatabaseParams()
    binning_input = assemblyStageChannel()
        .join(alignmentStageChannel())
        .map { sample, assembly_row, assembly_directory, alignment_row, alignment_directory ->
            tuple(sample, assembly_row, assembly_directory, alignment_directory)
        }
    BINNING(binning_input)
}

workflow MODULE7_REASSEMBLY {
    resolveDatabaseParams()
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

workflow MODULE8_SCAFFOLDING {
    resolveDatabaseParams()
    scaffolding_input = reassemblyStageChannel()
        .join(alignmentStageChannel())
        .join(libraryStageChannel('hic'))
        .map { sample, reassembly_row, reassembly_directory, alignment_row, alignment_directory, hic_row, hic_directory ->
            tuple(sample, reassembly_row, reassembly_directory, alignment_directory, hic_directory)
        }
    SCAFFOLDING(scaffolding_input)
}

workflow MODULE9_ANNOTATION {
    resolveDatabaseParams()
    ANNOTATION(reassemblyStageChannel())
}

workflow MODULE10_MGE {
    resolveDatabaseParams()

    if (params.mge_contact_dir) {
        mge_input = reassemblyStageChannel()
            .join(mgeContactStageChannel())
            .map { sample, reassembly_row, reassembly_directory, contact_row, contact_directory ->
                tuple(sample, reassembly_row, reassembly_directory, contact_directory)
            }
    } else {
        if (params.mge_alignment_dir) {
            mge_alignment_results = mgeAlignmentStageChannel()
        } else {
            mge_alignment_input = reassemblyStageChannel()
                .join(libraryStageChannel('hic'))
                .map { sample, reassembly_row, reassembly_directory, hic_row, hic_directory ->
                    tuple(sample, reassembly_row, reassembly_directory, hic_directory)
                }
            MGE_ALIGNMENT(mge_alignment_input)
            mge_alignment_results = MGE_ALIGNMENT.out.results
        }

        mge_contact_input = reassemblyStageChannel()
            .join(mge_alignment_results)
            .map { sample, reassembly_row, reassembly_directory, alignment_row, alignment_directory ->
                tuple(sample, reassembly_row, reassembly_directory, alignment_directory)
            }
        MGE_CONTACT(mge_contact_input)

        mge_input = reassemblyStageChannel()
            .join(MGE_CONTACT.out.results)
            .map { sample, reassembly_row, reassembly_directory, contact_row, contact_directory ->
                tuple(sample, reassembly_row, reassembly_directory, contact_directory)
            }
    }
    MGE(mge_input)
}

workflow {
    def selected = (params.entry_module ?: 'all').toString().trim().toLowerCase()

    if (selected in ['all', 'full', 'run_all']) {
        RUN_ALL()
    } else if (selected in ['module1', 'module1_preprocessing', 'preprocessing']) {
        MODULE1_PREPROCESSING()
    } else if (selected in ['module2', 'module2_assembly', 'assembly']) {
        MODULE2_ASSEMBLY()
    } else if (selected in ['module3', 'module3_alignment', 'alignment']) {
        MODULE3_ALIGNMENT()
    } else if (selected in ['module4', 'module4_coverage', 'coverage']) {
        MODULE4_COVERAGE()
    } else if (selected in ['module5', 'module5_contact', 'contact']) {
        MODULE5_CONTACT()
    } else if (selected in ['module6', 'module6_binning', 'binning']) {
        MODULE6_BINNING()
    } else if (selected in ['module7', 'module7_reassembly', 'reassembly']) {
        MODULE7_REASSEMBLY()
    } else if (selected in ['module8', 'module8_scaffolding', 'scaffolding']) {
        MODULE8_SCAFFOLDING()
    } else if (selected in ['module9', 'module9_annotation', 'annotation']) {
        MODULE9_ANNOTATION()
    } else if (selected in ['module10', 'module10_mge', 'mge']) {
        MODULE10_MGE()
    } else {
        error "Unknown --entry_module '${params.entry_module}'. Use all, module1, module2, module3, module4, module5, module6, module7, module8, module9, or module10."
    }
}
