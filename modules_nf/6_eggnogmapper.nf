#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process EGGNOGMAPPER {
    tag "Running eggNOG-mapper"
    label 'process_medium' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/08_eggnog-mapper", mode: 'copy', saveAs: { filename -> filename.replace("eggnog_results/", "") }

    input:
    path pangenome
    val db_path
    val type
    val mode
    val orthologs
    val go
    val pfam
    val diamond_sens

    output:
    path "eggnog_results/*", emit: eggnog_results

    script:
    """
    export EGGNOG_DATA_DIR=${db_path}

    awk '/^>/ {print ">" \$2; next} {print}' ${pangenome} > pan_genome_reference_pagoo.fa

    mkdir -p eggnog_results

    emapper.py \
        -i pan_genome_reference_pagoo.fa \
        --itype ${type} \
        -m ${mode} \
        --tax_scope Bacteria \
        --target_orthologs ${orthologs} \
        --go_evidence ${go} \
        --pfam_realign ${pfam} \
        --report_orthologs \
        --sensmode ${diamond_sens} \
        --output pan_genome_eggnog \
        --output_dir eggnog_results \
        --cpu ${task.cpus}
    """
}