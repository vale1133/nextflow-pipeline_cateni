#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process RAXML_TREE {
    tag "Building phylogenetic tree with RAxML"

    // Publish relevant outputs
    publishDir "${params.outdir}/8_raxml", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

    input:
    path alignment
    val runs

    output:
    path "RAxML_bestTree.core_gene.tree", emit: tree
    path "RAxML_info.core_gene.tree"

    script:
    """
    raxmlHPC \
        -m GTRCAT \
        -p 12345 \
        -s ${alignment} \
        -n core_gene.tree \
        -N ${runs} \
        -T ${task.cpus}
    """
}