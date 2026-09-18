#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Parameters ──────────────────────────────────────────────────────────────
params.samples    = null          // CSV: sample_id,sr1,sr2,lr
params.platonDb   = null
params.propagate  = 2
params.outdir     = "results"
params.useSpades     = false
params.perComponent  = false

// chopper defaults (matching Snakemake workflow)
params.chopper_minqual  = 9
params.chopper_minlen   = 500
params.chopper_headcrop = 75
params.chopper_tailcrop = 75

// ── Unknown parameter guard ────────────────────────────────────────────────
def knownParams = [
    'samples', 'platonDb', 'propagate', 'outdir',
    'useSpades', 'perComponent',
    'chopper_minqual', 'chopper_minlen', 'chopper_headcrop', 'chopper_tailcrop'
] as Set

// Nextflow injects these into params — skip them
def nextflowInternalParams = [
    'keep_work', 'help', 'version'
] as Set

def toKebab = { s -> s.replaceAll(/([A-Z])/, '-$1').toLowerCase() }

params.keySet().each { key ->
    if (!knownParams.contains(key) && !nextflowInternalParams.contains(key)) {
        def known = knownParams.sort().collect { "--${toKebab(it)}" }.join(', ')
        error "Unrecognized parameter: --${toKebab(key)}. Known parameters: ${known}"
    }
}

// ── Samplesheet validation ─────────────────────────────────────────────────

def validateSamplesheet(row) {
    def required = ['sample_id', 'sr1', 'sr2', 'lr']
    def missing  = required.findAll { !row.containsKey(it) || !row[it] }
    if (missing) {
        error "Samplesheet row missing required columns: ${missing.join(', ')}. Row: ${row}"
    }

    def files = ['sr1', 'sr2', 'lr']
    files.each { col ->
        def f = file(row[col])
        if (!f.exists()) {
            error "Sample '${row.sample_id}': file not found for '${col}': ${row[col]}"
        }
    }

    if (!row.sample_id.matches(/^[A-Za-z0-9._-]+$/)) {
        error "Sample '${row.sample_id}': sample_id contains invalid characters (use alphanumeric, dash, underscore, dot)"
    }

    return tuple(row.sample_id, file(row.sr1), file(row.sr2), file(row.lr))
}

// ── Processes ───────────────────────────────────────────────────────────────

process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(sr1), path(sr2), path(lr)

    output:
    tuple val(sample_id), path("${sample_id}.sr1.fastq.gz"),
                          path("${sample_id}.sr2.fastq.gz"), path(lr)
    path("${sample_id}.fastp.json"), emit: json

    script:
    """
    fastp --in1 ${sr1} --in2 ${sr2} \
          --out1 ${sample_id}.sr1.fastq.gz \
          --out2 ${sample_id}.sr2.fastq.gz \
          --unpaired1 ${sample_id}.unpaired.fastq.gz \
          --unpaired2 ${sample_id}.unpaired.fastq.gz \
          --json ${sample_id}.fastp.json \
          --thread ${task.cpus}
    """
}


process CHOPPER {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(sr1), path(sr2), path(lr)

    output:
    tuple val(sample_id), path(sr1), path(sr2), path("trim.lr.fastq.gz")

    script:
    """
    chopper \
        -q ${params.chopper_minqual} \
        --threads ${task.cpus} \
        -l ${params.chopper_minlen} \
        --headcrop ${params.chopper_headcrop} \
        --tailcrop ${params.chopper_tailcrop} \
        --input ${lr} \
        | pigz > trim.lr.fastq.gz
    """
}


process HYPLAS {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sr1), path(sr2), path(lr)
    path platon_db

    output:
    tuple val(sample_id), path("hyplas_out/*")

    script:
    def spades_flag     = params.useSpades    ? '--use-spades'     : ''
    def percomp_flag    = params.perComponent ? '--per-component'  : ''
    """
    hyplas \
        --platon-db ${platon_db} \
        -s ${sr1} ${sr2} \
        -l ${lr} \
        -o hyplas_out \
        -t ${task.cpus} \
        -p ${params.propagate} \
        --soft-fail \
        --keep-temp  \
        ${spades_flag} \
        ${percomp_flag} 
    """
}

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('reports/*')

    output:
    path("multiqc_report.html")
    path("multiqc_data")

    script:
    """
    multiqc --force reports/
    """
}

// ── Workflow ────────────────────────────────────────────────────────────────

workflow {
    if (!params.samples) {
        error "Please provide --samples <samples.csv>"
    }
    if (!params.platonDb) {
        error "Please provide --platon-db <path>"
    }

    Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row -> validateSamplesheet(row) }
        .set { samples_ch }

    // QC: short reads
    FASTP(samples_ch)

    // Trim long reads
    CHOPPER(FASTP.out[0])

    // Assembly
    platon_db_ch = file(params.platonDb, type: 'dir', checkIfExists: true)
    HYPLAS(CHOPPER.out, platon_db_ch)

    // Aggregate preprocessing QC reports
    FASTP.out.json
        .collect()
        | MULTIQC
}
