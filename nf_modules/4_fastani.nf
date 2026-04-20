#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FASTANI_ALL_VS_ALL {
    tag "ANI calculation with FastANI"
    label 'process_medium'  // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/4_fastani", pattern: "fastani_output", mode: 'copy'
    publishDir "${params.outdir}/4_fastani", pattern: "fastani_output.matrix", mode: 'copy'

    input:
    path genomes
    val files_extension

    output:
    path "fastani_output", emit: ani_results
    path "fastani_output.matrix", emit: ani_matrix, optional: true

    script:
    """
    find . -name "*.${files_extension}" > genomes_list.txt

    fastANI \
        --ql genomes_list.txt \
        --rl genomes_list.txt \
        -o fastani_output \
        -t ${task.cpus} \
        --matrix
    """
}