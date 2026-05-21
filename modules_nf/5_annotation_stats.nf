#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BAKTA_STATS {
    tag "Collecting Bakta annotation stats"

    // Publish relevant outputs
    publishDir "${params.outdir}/05_bakta/stats", mode: 'copy'

    input:
    path json_files
    path script_python

    output:
    path "annotation-stats.tsv", emit: stats

    script:
    """
    python ${script_python} *.json > annotation_stats.txt
    """
}