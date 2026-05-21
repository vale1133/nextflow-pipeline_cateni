#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CHECKM_QC {
    tag "Genome quality filtering with CheckM"
    label 'process_medium'  // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/03_checkm/filtered_genomes", pattern: "high_quality_genomes/*", mode: 'copy', saveAs: { filename -> filename.split('/').last() }
    publishDir "${params.outdir}/03_checkm", pattern: "quality_report.tsv", mode: 'copy'

    input:
    path genomes
    path db_path
    val files_extension
    val min_comp
    val max_cont

    output:
    path "quality_report.tsv", emit: quality_table
    path "high_quality_genomes/*", emit: filtered_genomes_checkm
    val count_genomes, emit: count

    script:
    count_genomes = genomes.size()
    """
    export CHECKM_DATA_PATH=${db_path}
    checkm lineage_wf \
        . \
        -t ${task.cpus} \
        -x ${files_extension} \
        -f quality_report.tsv \
        -r \
        checkm_output

    mkdir -p high_quality_genomes
    sed '/---/d' quality_report.tsv | sed 1d | awk '\$13>=${min_comp} && \$14<${max_cont} {print \$1}' | while read -r id; do
        file_match=\$(find . -maxdepth 1 -name "\${id}.${files_extension}" | head -n 1)
        if [ -n "\${file_match}" ]; then
            path=\$(readlink -f "\${file_match}")
            ln -s "\${path}" high_quality_genomes/
        fi
    done
    """
}