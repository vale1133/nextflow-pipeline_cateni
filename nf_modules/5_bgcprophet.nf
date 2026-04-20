#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BGC_PROPHET {
    tag "Running BGC-Prophet"

    maxForks 1

    // Publish relevant outputs
    publishDir "${params.outdir}/12_bgcprophet", mode: 'copy'

    input:
    path genomes
    path model_path
    path classifier_path
    val threshold
    val max_gap
    val min_count
    val classify_t
    val prefix

    output:
    path "*", emit: prophet_results

    script:
    """
    bgc_prophet pipeline --genomesDir . \
        --name ${prefix} \
        --modelPath ${model_path} \
        --threshold ${threshold} \
        --max_gap ${max_gap} \
        --min_count ${min_count} \
        --classifierPath ${classifier_path} \
        --classify_t ${classify_t}
    """
}