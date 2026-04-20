#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ROARY_PANGENOME_SPECIES {
    tag "Species pangenome calculation with Roary"
    label 'process_high' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/7.2_roary_species/${meta.species}", mode: 'copy'

    input:
    val meta
    path gffs
    val identity_sp

    output:
    path "gene_presence_absence.csv", emit: matrix
    path "summary_statistics.txt", emit: stats
    path "core_gene_alignment.aln", emit: alignment
    path "*.Rtab", emit: rtabs
    path "Rplots.pdf", emit: rplots
    path "pan_genome_reference.fa", emit: pangenome_ref
    path "*.png", emit: pangenome_plots

    script:
    """
    roary \
        -p ${task.cpus} \
        -e -n \
        -v \
        -r \
        -i ${identity_sp} \
        *.gff
    """
}