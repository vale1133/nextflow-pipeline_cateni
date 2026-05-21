#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BGC_PROPHET {
    tag "Running BGC-Prophet"

    scratch true

    // Publish relevant outputs
    publishDir "${params.outdir}/11_bgcprophet", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

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
    path "prophet_output/${prefix}*", emit: prophet_results

    script:
    """
    bgc_prophet pipeline --genomesDir . \
        --name ${prefix} \
        --modelPath ${model_path} \
        --threshold ${threshold} \
        --max_gap ${max_gap} \
        --min_count ${min_count} \
        --classifierPath ${classifier_path} \
        --classify_t ${classify_t} \
        --outputPath prophet_output
    """
}