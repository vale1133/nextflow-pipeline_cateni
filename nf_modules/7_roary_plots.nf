#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ROARY_PLOTS {
    tag "Generating Roary plots"
    label 'process_low' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/7_roary/roary_plots", mode: 'copy'

    input:
    path tree_file
    path matrix_file
    path script_python

    output:
    path "*", emit: plots

    script:
    """
    python ${script_python} ${tree_file} ${matrix_file}
    """
}