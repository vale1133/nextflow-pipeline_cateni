#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process RENAME_FAA {
    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotator}.faa")

    script:
    """
    ln -s $faa ${meta.id}_${meta.annotator}.faa
    """
}
