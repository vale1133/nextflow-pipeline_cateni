#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { GTDB_FILTER } from './modules_nf/1_gtdbtk.nf'
include { DEREPLICATE_SEQS } from './modules_nf/2_rmdup.nf'
include { CHECKM_QC } from './modules_nf/3_checkm.nf'
include { FASTANI_ALL_VS_ALL } from './modules_nf/4_fastani.nf'
include { BAKTA_ANNOTATION } from './modules_nf/4_bakta.nf'
include { BAKTA_STATS } from './modules_nf/4-1_annotation_stats.nf'
include { PROKKA_ANNOTATION } from './modules_nf/4_prokka.nf'
include { ROARY_PANGENOME } from './modules_nf/5_roary.nf'
include { RAXML_TREE } from './modules_nf/6_raxml.nf'
include { ROARY_PLOTS } from './modules_nf/7_roary_plots.nf'
include { EGGNOGMAPPER } from './modules_nf/6_eggnogmapper.nf'
include { DEEPBGC } from './modules_nf/4_deepbgc.nf'
include { ANTISMASH } from './modules_nf/4_antismash.nf'
include { BGC_PROPHET } from './modules_nf/5_bgcprophet.nf'
include { RENAME_FAA } from './modules_nf/rename_faa.nf'
include { ROARY_PANGENOME_SPECIES } from './modules_nf/5_roary_especies.nf'

workflow {
    // 1. DEFINE CHANNELS
    // 1.1. GTDB-Tk database path channel
    gtdbtk_db_ch = channel.value(params.gtdbtk_db)

    // 1.2. Genomes path channel
    genomes_ch = channel.fromPath(params.genomes, type: 'dir')

    // 1.3. Files extension channel
    extension_ch = channel.value(params.files_extension)

    // 1.4. CheckM database path channel
    checkm_db_ch = channel.value(params.checkm_db)

    // 1.5. Genome completeness and contamination channels
    comp_ch = channel.value(params.checkm_min_comp)
    cont_ch = channel.value(params.checkm_max_cont)

    // 1.6. Taxonomy channel
    tax_ch = channel.value(params.tax)

    // 1.7. Bakta database path channel
    bakta_db_ch = channel.value(params.bakta_db)

    // 1.8. Bakta stats script path channel
    bakta_script_ch = channel.fromPath(params.bakta_stats_script)

    // 1.9. Roary minimum identity channel
    roary_identity_ch = channel.value(params.roary_min_identity)

    // 1.10. RAxML runs channel
    raxml_runs_ch = channel.value(params.raxml_runs)

    // 1.11. Roary plots script path channel
    roary_plots_ch = channel.fromPath(params.roary_plots_script)

    // 1.12. EggNOG-mapper database path channel
    eggnog_db_ch = channel.value(params.eggnog_db)

    // 1.13. EggNOG-mapper input type channel
    eggnog_input_type_ch = channel.value(params.eggnog_input_type)

    // 1.14. EggNOG-mapper search mode channel
    eggnog_search_mode_ch = channel.value(params.eggnog_search_mode)

    // 1.15. EggNOG-mapper orthologs channel
    eggnog_orthologs_ch = channel.value(params.eggnog_orthologs)

    // 1.16. EggNOG-mapper GO evidence channel
    eggnog_go_evidence_ch = channel.value(params.eggnog_go_evidence)

    // 1.17. EggNOG-mapper PFAM realign channel
    eggnog_pfam_realign_ch = channel.value(params.eggnog_pfam_realign)

    // 1.18. EggNOG-mapper sensitivity mode channel
    eggnog_sens_mode_ch = channel.value(params.eggnog_sens_mode)

    // 1.19. DeepBGC models path channel
    deepbgc_models_ch = channel.value(params.deepbgc_models)

    // 1.20. DeepBGC model channel
    deepbgc_model_ch = channel.value(params.deepbgc_model)

    // 1.21. DeepBGC extracting score channel
    deepbgc_extracting_score_ch = channel.value(params.deepbgc_extracting_score)

    // 1.22. DeepBGC classification score channel
    deepbgc_classification_score_ch = channel.value(params.deepbgc_classification_score)

    // 1.23. DeepBGC protein gap channel
    deepbgc_prot_gap_ch = channel.value(params.deepbgc_prot_gap)

    // 1.24. DeepBGC nucleotide gap channel
    deepbgc_nucl_gap_ch = channel.value(params.deepbgc_nucl_gap)

    // 1.25. antiSMASH gene finder tool channel
    antismash_genefinder_ch = channel.value(params.antismash_genefinder)

    // 1.26. BGC-Prophet model path channel
    prophet_model_ch = channel.value(params.prophet_model)

    // 1.27. BGC-Prophet classifier path channel
    prophet_classifier_ch = channel.value(params.prophet_classifier)

    // 1.28. BGC-Prophet prediction threshold channel
    prophet_pred_t_ch = channel.value(params.prophet_pred_t)

    // 1.29. BGC-Prophet gene gap threshold channel
    prophet_gene_gap_ch = channel.value(params.prophet_gene_gap)

    // 1.30. BGC-Prophet minimum number of genes channel
    prophet_min_genes_ch = channel.value(params.prophet_min_genes)

    // 1.31. BGC-Prophet classification threshold channel
    prophet_class_t_ch = channel.value(params.prophet_class_t)

    // 1.32. BGC-Prophet name channel
    prophet_prefix_ch = channel.value(params.prophet_prefix)

    // 1.33. Input channel for Bakta and Prokka annotation processes, and for DeepBGC and antiSMASH
    input_ch = CHECKM_QC.out.filtered_genomes_checkm
        .flatten()
        .map { file -> 
            def sample_id = file.name.minus(".${params.files_extension}")
            def meta = [ id: sample_id ] 
            return [ meta, file ]
        }

    // 1.34. antiSMASH database path channel
    antismash_db_ch = channel.value(params.antismash_db)

    // 1.35. Define Roary for species input channel
    species_ch = channel
        .fromPath("${params.outdir}/1_gtdbtk/gtdbtk.summary.tsv")
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def tax = row.classification.split(';')
            def species = tax.find { it.startsWith('s__') } ?: 'Unknown species'
            return [ row.user_genome, species.replace('s__', '').replace(' ', '_') ]
        }

    ch_prokka_pre_join = PROKKA_ANNOTATION.out.gff.map { meta, gff -> [ meta.id, gff ] }

    ch_gffs_w_species = ch_prokka_pre_join
        .join(species_ch)
        .map { id, gff, species ->
            def new_meta = [id: id, species: species]
            return [ new_meta, gff ]
        }

    ch_branched_gffs = ch_gffs_w_species.branch { meta, gff ->
        mitsuokai: meta.species == 'Catenibacterium_mitsuokai'
        mitsuokai_a: meta.species == 'Catenibacterium_mitsuokai_A'
        }

    // 1.36. Roary minimum identity channel for species pangenome calculation
    identity_sp_ch = channel.value(params.roary_min_identity_sp)


    // 2. RUN GTDBTK
    GTDB_FILTER( gtdbtk_db_ch, genomes_ch, extension_ch, tax_ch )

    // 2.1. Count genomes before and after filtering
    def before_count_gtdbtk = genomes_ch
        .map { dir -> 
            dir.listFiles { file -> file.name.endsWith(".${params.files_extension}") }
            .size()
        }

    def after_count_gtdbtk = GTDB_FILTER.out.filtered_genomes_gtdbtk
        .flatten()
        .count()

    // 2.2. Print results
    GTDB_FILTER.out.summary
        .combine(before_count_gtdbtk)
        .combine(after_count_gtdbtk)
        .subscribe { summary, before, after ->
            log.info "------------------------------------------------------------------------"
    // 2.2.1. Print results path and directory
            log.info "📄 GTDB-TK Summary in ${params.outdir}/1_gtdbtk/gtdbtk_summary.tsv"
            log.info "📂 GTDB-TK Filtered genomes in ${params.outdir}/1_gtdbtk/filtered_genomes"
    // 2.2.2. See how many genomes passed the filtering step
            log.info "🧬 GTDB-TK Genomes before filtering step: ${before}"
            log.info "🧬 GTDB-TK Genomes after filtering step: ${after} (${params.tax})"
            log.info "------------------------------------------------------------------------"
        }


    // 3. RUN SEQKIT RMDUP
    DEREPLICATE_SEQS( GTDB_FILTER.out.filtered_genomes_gtdbtk.collect() )

    // 3.1. Count genomes before and after dereplication
    def before_count_seqkit = DEREPLICATE_SEQS.out.count

    def after_count_seqkit = DEREPLICATE_SEQS.out.unique_genomes
        .flatten()
        .count()

    // 3.2. Print results
    DEREPLICATE_SEQS.out.report
        .combine(before_count_seqkit)
        .combine(after_count_seqkit)
        .subscribe { report, before, after ->
            log.info "------------------------------------------------------------------------"
    // 3.2.1. Print results path and directory
            log.info "📄 SEQKIT Report in ${params.outdir}/2_seqkit/rmdup_report.txt"
            log.info "📄 SEQKIT Duplicated sequences in ${params.outdir}/2_seqkit/dup_sequences.fasta"
            log.info "📂 SEQKIT Filtered genomes in ${params.outdir}/2_seqkit/filtered_genomes"
    // 3.2.2. See how many genomes passed the filtering step
            log.info "🧬 SEQKIT Genomes before duplicate removal: ${before}"
            log.info "🧬 SEQKIT Genomes after duplicate removal: ${after}"
            log.info "------------------------------------------------------------------------"
        }


    // 4. RUN CHECKM
    CHECKM_QC( DEREPLICATE_SEQS.out.unique_genomes.collect(), checkm_db_ch, extension_ch, comp_ch, cont_ch )
    
    // 4.1. Count genomes before and after dereplication
    def before_count_checkm = CHECKM_QC.out.count

    def after_count_checkm = CHECKM_QC.out.filtered_genomes_checkm
        .flatten()
        .count()

    // 4.2. Print results
    CHECKM_QC.out.quality_table
        .combine(before_count_checkm)
        .combine(after_count_checkm)
        .subscribe { table, before, after ->
            log.info "------------------------------------------------------------------------"
    // 4.2.1. Print results path and directory
            log.info "📄 CHECKM Table in ${params.outdir}/3_checkm/checkm_results.tsv"
            log.info "📂 CHECKM Filtered genomes in ${params.outdir}/3_checkm/filtered_genomes"
    // 4.2.2. See how many genomes passed the filtering step
            log.info "🧬 CHECKM Genomes before filtering step: ${before}"
            log.info "🧬 CHECKM Genomes after filtering step: ${after}"
            log.info "------------------------------------------------------------------------"
    }


    // 5. RUN FASTANI
    FASTANI_ALL_VS_ALL( CHECKM_QC.out.filtered_genomes_checkm.collect(), extension_ch )

    // 5.1. Print results paths
    FASTANI_ALL_VS_ALL.out.ani_results
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 FASTANI Results in ${params.outdir}/4_fastani"
            log.info "------------------------------------------------------------------------"
        }


    // 6. RUN BAKTA 
    BAKTA_ANNOTATION( input_ch, bakta_db_ch, tax_ch )

    // 6.2. Print Bakta results
    BAKTA_ANNOTATION.out.gff
        .collect()
        .subscribe { gffs ->
            log.info "------------------------------------------------------------------------"
    // 6.2.1. Print results directory
            log.info "📂 BAKTA Results in ${params.outdir}/5_bakta"
    // 6.2.2. See how many genomes were annotated
            log.info "🧬 BAKTA ${gffs.size()} genomes were correctly annotated by Bakta"
            log.info "------------------------------------------------------------------------"
    }

    // 6.4. RUN COLLECT ANNOTATION STATS
    BAKTA_STATS( BAKTA_ANNOTATION.out.json.collect(), bakta_script_ch )

    // 6.5. Print Bakta stats results directory
    BAKTA_STATS.out.stats
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 BAKTA ANNOTATION STATS Results in ${params.outdir}/5_bakta/stats"
            log.info "------------------------------------------------------------------------"
        }


    // 7. RUN PROKKA
    PROKKA_ANNOTATION( input_ch, tax_ch )

    // 7.2. Print results
    PROKKA_ANNOTATION.out.gff
        .collect()
        .subscribe { gffs ->
            log.info "------------------------------------------------------------------------"
    // 7.2.1. Print results directory
            log.info "📂 PROKKA Results in ${params.outdir}/6_prokka"
    // 7.2.2. See how many genomes were annotated
            log.info "🧬 PROKKA ${gffs.size()} genomes were correctly annotated by Prokka"
            log.info "------------------------------------------------------------------------"
    }


    // 8. RUN ROARY
    ROARY_PANGENOME( PROKKA_ANNOTATION.out.gff.collect(), roary_identity_ch )

    // 8.1. Print results directory
    ROARY_PANGENOME.out.stats
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 ROARY Results in ${params.outdir}/7_roary"
            log.info "------------------------------------------------------------------------"
        }

    // 8.2. RUN ROARY FOR SPECIES
    ROARY_PANGENOME_SPECIES(
        ch_branched_gffs.mitsuokai.map{ meta, gff -> meta.species }.collect().first(), 
        ch_branched_gffs.mitsuokai.map{ meta, gff -> gff }.collect(),
        identity_sp_ch
    )

    ROARY_PANGENOME_SPECIES(
        ch_branched_gffs.mitsuokai_a.map{ meta, gff -> meta.species }.collect().first(), 
        ch_branched_gffs.mitsuokai_a.map{ meta, gff -> gff }.collect(),
        identity_sp_ch
    )


    // 9. RUN RAXML
    RAXML_TREE( ROARY_PANGENOME.out.alignment, raxml_runs_ch )

    // 9.1. Print results directory
    RAXML_TREE.out.tree
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 RAXML Results in ${params.outdir}/8_raxml"
            log.info "------------------------------------------------------------------------"
        }


    // 10. RUN ROARY PLOTS
    ROARY_PLOTS( RAXML_TREE.out.tree, ROARY_PANGENOME.out.matrix, roary_plots_ch )
    
    // 10.1. Print results directory
    ROARY_PLOTS.out.plots
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 ROARY PLOTS Results in ${params.outdir}/7_roary/roary_plots"
            log.info "------------------------------------------------------------------------"
        }


    // 11. RUN EGGNOG-MAPPER
    EGGNOGMAPPER( ROARY_PANGENOME.out.pangenome_ref, eggnog_db_ch, eggnog_input_type_ch, eggnog_search_mode_ch, eggnog_orthologs_ch, eggnog_go_evidence_ch, eggnog_pfam_realign_ch, eggnog_sens_mode_ch )

    // 11.1. Print results directory
    EGGNOGMAPPER.out.eggnog_results
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 EGGNOG-MAPPER Results in ${params.outdir}/9_eggnog-mapper"
            log.info "------------------------------------------------------------------------"
        }


    // 12. RUN DEEPBGC
    DEEPBGC( input_ch, deepbgc_models_ch, deepbgc_model_ch, deepbgc_extracting_score_ch, deepbgc_classification_score_ch, deepbgc_prot_gap_ch, deepbgc_nucl_gap_ch )

    // 12.2. Print results directory
    DEEPBGC.out.deepbgc_results
        .last()
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 DEEPBGC Results in ${params.outdir}/10_deepbgc"
            log.info "------------------------------------------------------------------------"
        }


    // 13. RUN ANTISMASH
    ANTISMASH( input_ch, antismash_genefinder_ch, antismash_db_ch )

    // 13.2. Print results directory
    ANTISMASH.out.antismash_results
        .last()
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 ANTISMASH Results in ${params.outdir}/11_antismash"
            log.info "------------------------------------------------------------------------"
        }


    // 14. BGC-PROPHET
    // 14.1. Define BGC-Prophet input channel
    prokka_faa_ch = PROKKA_ANNOTATION.out.faa
        .map { meta, faa -> 
            [ meta + [annotator: 'prokka'], faa ] 
    }

    bakta_faa_ch = BAKTA_ANNOTATION.out.faa
        .branch { meta, faa ->
            hypotheticals: faa.name.endsWith('.hypotheticals.faa')
                return [ meta + [annotator: 'bakta_hypotheticals'], faa ]
            normal: faa.name.endsWith('.faa')
                return [ meta + [annotator: 'bakta'], faa ]
        }

    combined_faa_ch = prokka_faa_ch
        .mix(bakta_faa_ch.normal)
        .mix(bakta_faa_ch.hypotheticals)

    all_faa_ch = RENAME_FAA(combined_faa_ch)
        .map { meta, faa -> faa }
        .collect()

    // 14.2. RUN BGC-PROPHET
    BGC_PROPHET( all_faa_ch, prophet_model_ch, prophet_classifier_ch, prophet_pred_t_ch, prophet_gene_gap_ch, prophet_min_genes_ch, prophet_class_t_ch, prophet_prefix_ch )

    // 14.3. Print results directory
    BGC_PROPHET.out.prophet_results
        .subscribe { 
            log.info "------------------------------------------------------------------------"
            log.info "📂 BGC-PROPHET Results in ${params.outdir}/12_bgcprophet"
            log.info "------------------------------------------------------------------------"
        }
}