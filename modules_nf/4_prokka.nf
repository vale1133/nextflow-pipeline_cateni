#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PROKKA_ANNOTATION {
    tag "Prokka annotation: ${meta.id}"
    label 'process_medium' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/05_prokka/${meta.id}", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

    input:
    tuple val(meta), path(genome)
    val taxa

    output:
    tuple val(meta), path("prokka_output/*.gff"), emit: gff
    tuple val(meta), path("prokka_output/${meta.id}.faa"), emit: proteins
    path("prokka_output/*.ffn"), emit: genes
    path("prokka_output/*.txt"), emit: stats
    
    script:
    """
    prokka \
        --outdir prokka_output \
        --prefix ${meta.id} \
        --locustag ${meta.id} \
        --cpus ${task.cpus} \
        --kingdom Bacteria \
        --genus ${taxa} \
        ${genome}
    """
}