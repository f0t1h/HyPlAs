#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Parameters ──────────────────────────────────────────────────────────────
params.samples    = null          // CSV: sample_id,sr1,sr2,lr
params.platonDb   = null
params.propagate  = 0
params.outdir     = "results"
params.use_spades = false

// chopper defaults (matching Snakemake workflow)
params.chopper_minqual  = 9
params.chopper_minlen   = 500
params.chopper_headcrop = 75
params.chopper_tailcrop = 75

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

process NANOPLOT_RAW {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc/nanoplot_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(sr1), path(sr2), path(lr)

    output:
    path("${sample_id}_raw_NanoStats.txt"), emit: stats

    script:
    """
    NanoPlot --fastq ${lr} \
             --threads ${task.cpus} \
             --prefix ${sample_id}_raw_ \
             --no_static
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

process NANOPLOT_TRIMMED {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc/nanoplot_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(sr1), path(sr2), path(lr)

    output:
    path("${sample_id}_trimmed_NanoStats.txt"), emit: stats

    script:
    """
    NanoPlot --fastq ${lr} \
             --threads ${task.cpus} \
             --prefix ${sample_id}_trimmed_ \
             --no_static
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
    def spades_flag = params.use_spades ? '--use-spades' : ''
    """
    hyplas \
        --platon-db ${platon_db} \
        -s ${sr1} ${sr2} \
        -l ${lr} \
        -o hyplas_out \
        -t ${task.cpus} \
        -p ${params.propagate} \
        --soft-fail \
        ${spades_flag}
    """
}

process QUAST {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly_files)

    output:
    path("${sample_id}_quast"), emit: report

    script:
    """
    # Run QUAST on the final iteration assembly
    FINAL_FASTA=\$(ls -1 plasmids.final.it*.fasta 2>/dev/null | sort -t. -k4 -n | tail -1)

    if [ -n "\$FINAL_FASTA" ] && [ -s "\$FINAL_FASTA" ]; then
        quast \$FINAL_FASTA \
            --output-dir ${sample_id}_quast \
            --label ${sample_id} \
            --min-contig 0 \
            --threads ${task.cpus}
    else
        # Empty assembly — create minimal report so pipeline doesn't fail
        mkdir -p ${sample_id}_quast
        echo "No plasmid contigs assembled for ${sample_id}" > ${sample_id}_quast/report.txt
    fi
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

    // QC: long reads (raw) — runs in parallel with FASTP
    NANOPLOT_RAW(samples_ch)

    // Trim long reads
    CHOPPER(FASTP.out[0])

    // QC: long reads (trimmed)
    NANOPLOT_TRIMMED(CHOPPER.out)

    // Assembly
    platon_db_ch = file(params.platonDb, type: 'dir', checkIfExists: true)
    HYPLAS(CHOPPER.out, platon_db_ch)

    // Assembly QC
    QUAST(HYPLAS.out)

    // Aggregate all QC reports
    FASTP.out.json
        .mix(
            NANOPLOT_RAW.out.stats,
            NANOPLOT_TRIMMED.out.stats,
            QUAST.out.report
        )
        .collect()
        | MULTIQC
}
