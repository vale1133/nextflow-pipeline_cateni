#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ANTISMASH {
    tag "Running antiSMASH: ${meta.id}"

    // Publish relevant outputs
    publishDir "${params.outdir}/11_antiSMASH/${meta.id}", mode: 'copy', saveAs: { filename -> filename.replace("antismash/", "") }

    input:
    tuple val(meta), path(genome)
    val gene_tool
    path db_path

    output:
    path "antismash/*", emit: antismash_results

    script:
    """
    antismash \
        --output-dir antismash \
        --output-basename "antismash_${meta.id}" \
        -t bacteria \
        --cb-knownclusters \
        --cc-mibig \
        --asf \
        --genefinding-tool ${gene_tool} \
        --rre \
        --pfam2go \
        --fullhmmer \
        --databases ${db_path} \
        ${genome}
    """
}