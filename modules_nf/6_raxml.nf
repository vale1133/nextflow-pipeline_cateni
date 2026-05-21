#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process RAXML_TREE {
    tag "Building phylogenetic tree with RAxML"

    // Publish relevant outputs
    publishDir "${params.outdir}/07_raxml", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

    input:
    path alignment
    val runs

    output:
    path "RAxML_bestTree.core_gene.tree_bs", emit: tree
    path "*.core_gene.tree_bs"

    script:
    """
    raxmlHPC \
        -f a \
        -p 12345 \
        -s ${alignment} \
        -x 12345 \
        -N ${runs} \
        -m GTRCAT \
        -n core_gene.tree_bs \
        -T ${task.cpus}
    """
}