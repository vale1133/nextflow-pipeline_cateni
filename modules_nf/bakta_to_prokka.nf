#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BAKTA_TO_PROKKA_GFF {
    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("${meta.id}.gff")

    script:
    """
    cp "${gff3}" "${meta.id}.gff"
    """
}