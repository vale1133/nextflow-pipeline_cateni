#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DEREPLICATE_SEQS {
    tag "Removing duplicated genomes with SeqKit"
    
    // Publish relevant outputs
    publishDir "${params.outdir}/2_seqkit/filtered_genomes", pattern: "unique_genomes/*", mode: 'copy', saveAs: { filename -> filename.split('/').last() }
    publishDir "${params.outdir}/2_seqkit", pattern: "rmdup_report.txt", mode: 'copy'
    publishDir "${params.outdir}/2_seqkit", pattern: "dup_sequences.fasta", mode: 'copy'

    input:
    path genomes

    output:
    path "unique_genomes/*", emit: unique_genomes
    path "rmdup_report.txt", emit: report
    path "dup_sequences.fasta", emit: duplicates
    val count_genomes, emit: count

    script:
    count_genomes = genomes.size()
    """
    mkdir -p unique_genomes
    for fasta in ${genomes}; do
        filename=\$(basename \$fasta)
        echo ">\$filename" >> all_genomes_linear.fasta
        seqkit sort -n \$fasta | seqkit seq -w 0 | grep -v "^>" | tr -d "\n" >> all_genomes_linear.fasta
        echo "" >> all_genomes_linear.fasta
    done

    seqkit rmdup -s -i all_genomes_linear.fasta \
        -o clean_sequences.fasta \
        -D rmdup_report.txt \
        -d dup_sequences.fasta

    grep "^>" clean_sequences.fasta | sed 's/>//' > unique_genomes_list.txt

    while read -r file_to_save; do
        file_match=\$(find . -maxdepth 1 -name "\${file_to_save}" | head -n 1)
        if [ -n "\${file_match}" ]; then
            path=\$(readlink -f "\${file_match}")
            ln -s "\${path}" unique_genomes/
        fi
    done < unique_genomes_list.txt
    """
}