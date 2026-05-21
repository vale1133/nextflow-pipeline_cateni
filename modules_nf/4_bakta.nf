#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BAKTA_ANNOTATION {
    tag "Bakta annotation: ${meta.id}"
    label 'process_medium' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/05_bakta/${meta.id}", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

    input:
    tuple val(meta), path(genome)
    path db_path
    val taxa

    output:
    tuple val(meta), path("bakta_output/*.gff3"), emit: gff
    path "bakta_output/*.json", emit: json
    tuple val(meta), path("bakta_output/*.faa"), emit: proteins
    path "bakta_output/*.ffn", emit: genes
    path "bakta_output/*.gbff", emit: genbank
    path "bakta_output/*.txt", emit: summary
    path "bakta_output/*.tsv", emit: tsv

    script:
    """
    bakta \
        --db ${db_path} \
        --output bakta_output \
        --prefix ${meta.id} \
        --threads ${task.cpus} \
        --genus ${taxa} \
        ${genome}
    """
}