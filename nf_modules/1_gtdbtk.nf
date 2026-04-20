#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GTDB_FILTER {
    tag "Filtering ${params.tax} genomes with GTDB-Tk"
    label 'process_high' // CPU/RAM configuration in config file

    // Publish relevant outputs
    publishDir "${params.outdir}/1_gtdbtk", pattern: "gtdbtk_summary.tsv", mode: 'copy'
    publishDir "${params.outdir}/1_gtdbtk/filtered_genomes", pattern: "selected/*", mode: 'copy', saveAs: { filename -> filename.split('/').last() }

    input:
    path db_path
    path genomes_dir
    val files_extension
    val taxa

    output:
    path "selected/*", emit: filtered_genomes_gtdbtk
    path "gtdbtk_summary.tsv", emit: summary

    script:
    """
    export GTDBTK_DATA_PATH=${db_path}

    gtdbtk classify_wf \\
        --genome_dir ${genomes_dir} \\
        --out_dir gtdbtk_output \\
        --extension ${files_extension} \\
        --cpus ${task.cpus}

    mv gtdbtk_output/classify/gtdbtk.bac120.summary.tsv gtdbtk_summary.tsv

    mkdir -p selected
    sed 1d gtdbtk_summary.tsv | awk -F '\\t' -v pat="${taxa}" 'tolower(\$2) ~ tolower(pat) {print \$1}' | while read -r id; do
        file_match=\$(find "${genomes_dir}"/* -maxdepth 1 -name "\${id}*" | head -n 1)
        if [ -n "\${file_match}" ]; then
            path=\$(readlink -f "\${file_match}")
            ln -s "\${path}" selected/
        fi
    done
    """
}