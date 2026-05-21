#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ROARY_PANGENOME_SPECIES {
    tag "Species pangenome calculation with Roary: ${species}"
    label 'process_high' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/06_roarySpecies/${species}", mode: 'copy'

    input:
    tuple val(species), path(gffs)
    val identity_sp

    output:
    path "gene_presence_absence.csv", emit: matrix
    path "summary_statistics.txt", emit: stats
    path "core_gene_alignment.aln", emit: alignment
    path "*.Rtab", emit: rtabs
    path "Rplots.pdf", emit: rplots
    path "pan_genome_reference.fa", emit: pangenome_ref

    script:
    """
    roary \
        -p ${task.cpus} \
        -e -n \
        -r \
        -i ${identity_sp} \
        *.gff
    """
}