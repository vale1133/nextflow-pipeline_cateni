#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CHECKM_QC {
    tag "Genome quality filtering with CheckM"
    label 'process_medium'  // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/3_checkm/filtered_genomes", pattern: "high_quality_genomes/*", mode: 'copy', saveAs: { filename -> filename.split('/').last() }
    publishDir "${params.outdir}/3_checkm", pattern: "checkm_results.tsv",mode: 'copy'

    input:
    path genomes
    path db_path
    val files_extension
    val min_comp
    val max_cont

    output:
    path "checkm_results.tsv", emit: quality_table
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
        -f checkm_results.tsv \
        -r \
        checkm_output

    mkdir -p high_quality_genomes
    sed '/---/d' checkm_results.tsv | sed 1d | awk '\$13>=${min_comp} && \$14<${max_cont} {print \$1}' | while read -r id; do
        file_match=\$(find . -maxdepth 1 -name "\${id}.${files_extension}" | head -n 1)
        if [ -n "\${file_match}" ]; then
            path=\$(readlink -f "\${file_match}")
            ln -s "\${path}" high_quality_genomes/
        fi
    done
    """
}