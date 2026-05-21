#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DEEPBGC {
    tag "Running DeepBGC: ${meta.id}"
    label 'process_medium' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/09_deepbgc/${meta.id}", mode: 'copy', saveAs: { filename -> filename.replace("deepbgc_output/", "") }

    input:
    tuple val(meta), path(genome)
    path models_path
    val model
    val extracting_score
    val classification_score
    val prot_gap
    val nucl_gap

    output:
    path "deepbgc_output/deepbgc_output*", emit: deepbgc_results
    path "deepbgc_output/evaluation/deepbgc_output*", emit: images

    script:
    """
    export DEEPBGC_DOWNLOADS_DIR=${models_path}

    deepbgc pipeline ${genome} \
        --output deepbgc_output \
        --prodigal-meta-mode \
        -d ${model} \
        -s ${extracting_score} \
        --merge-max-protein-gap ${prot_gap} \
        --merge-max-nucl-gap ${nucl_gap} \
        -c product_class -c product_activity \
        --classifier-score ${classification_score}
    """
}