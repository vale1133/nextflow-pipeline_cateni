#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CHECKM2_QC {
    tag "Genome quality filtering with CheckM2"
    label 'process_high'  // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/03_checkm2/filtered_genomes", pattern: "high_quality_genomes/*", mode: 'copy', saveAs: { filename -> filename.split('/').last() }
    publishDir "${params.outdir}/03_checkm2", pattern: "checkm2_output/quality_report.tsv", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

    input:
    path genomes
    path db_path
    val files_extension
    val min_comp
    val max_cont

    output:
    path "checkm2_output/quality_report.tsv", emit: quality_table
    path "high_quality_genomes/*", emit: filtered_genomes_checkm
    val count_genomes, emit: count

    script:
    count_genomes = genomes.size()
    """
    export CHECKM2DB=${db_path}
    checkm2 predict \
        --threads ${task.cpus} \
        --input . \
        --output-directory checkm2_output

    mkdir -p high_quality_genomes
    sed '/---/d' checkm2_output/quality_report.tsv | sed 1d | awk '\$2>=${min_comp} && \$3<${max_cont} {print \$1}' | while read -r id; do
        file_match=\$(find . -maxdepth 1 -name "\${id}.${files_extension}" | head -n 1)
        if [ -n "\${file_match}" ]; then
            path=\$(readlink -f "\${file_match}")
            ln -s "\${path}" high_quality_genomes/
        fi
    done
    """
}