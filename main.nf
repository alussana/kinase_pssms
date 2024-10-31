#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { dl_human } from './modules/uniprot_id_dict'
include { uniprotac_to_ensp_dict } from './modules/uniprot_id_dict'
include { IDa2uniprot2IDb } from './modules/uniprot_id_dict'

include { dl_string_aliases } from './modules/string_id_dict'

include { get_ser_thr_kinome_2023_suppl_table_2 } from './modules/pssm'
include { get_ser_thr_kinase_pssm } from './modules/pssm'
include { get_tyr_kinome_2024_suppl_table_2 } from './modules/pssm'
include { get_tyr_kinase_pssm } from './modules/pssm'
include { s_t_pssm_h5_alias_mapping } from './modules/pssm'
include { y_pssm_h5_alias_mapping } from './modules/pssm'

include { psp_kinase_substrates } from './modules/phosphositeplus'
include { kin_phos_clusters } from './modules/phosphositeplus'
include { translate_kin_phos_clusters } from './modules/phosphositeplus'
include { hs_phosphoproteome } from './modules/phosphositeplus'

include { s_t_pssm_background_scores } from './modules/bg_scores'
include { y_pssm_background_scores } from './modules/bg_scores'

include { publish } from './modules/utils'
include { split } from './modules/utils'
include { concatenate } from './modules/utils'


workflow GET_STRING_ID_DICT {

    main:
        string_id_dict = dl_string_aliases()

    emit:
        string_id_dict

}


workflow SER_THR_KINASES_PSSM {

    take:
        string_id_dict

    main:
        kinome_2023_suppl_table_2 = get_ser_thr_kinome_2023_suppl_table_2()
        ser_thr_pssms = get_ser_thr_kinase_pssm( kinome_2023_suppl_table_2 )
        translated_ser_thr_pssms = s_t_pssm_h5_alias_mapping( ser_thr_pssms.pssm_dict_h5,
                                                          string_id_dict ).pssm_dict_h5

    emit:
        translated_ser_thr_pssms

}


workflow TYR_KINASES_PSSM {

    take:
        string_id_dict

    main:
        kinome_2024_suppl_table_2 = get_tyr_kinome_2024_suppl_table_2()
        tyr_pssms = get_tyr_kinase_pssm( kinome_2024_suppl_table_2)
        translated_tyr_pssms = y_pssm_h5_alias_mapping( tyr_pssms.pssm_dict_h5,
                                                      string_id_dict ).pssm_dict_h5

    emit:
        translated_tyr_pssms

}


workflow PSP {

    take:
        psp_kinase_substrates_dataset_gz
        psp_phosphosites_dataset_gz
        string_id_dict

    main:
        kinase_substrates = psp_kinase_substrates( psp_kinase_substrates_dataset_gz )
        kin_sub_clusters_untranslated = kin_phos_clusters( kinase_substrates )
        kin_sub_clusters = translate_kin_phos_clusters( kin_sub_clusters_untranslated,
                                                        string_id_dict )
        human_phosphosites = hs_phosphoproteome( psp_phosphosites_dataset_gz )

    emit:
        kin_sub_clusters
        human_phosphosites

}


workflow S_T_PSSM_BACKGROUND_SCORES {

    take:
        human_phosphosites_tsv
        pssm_dict_h5

    main:
        pssm_bg_scores = s_t_pssm_background_scores( human_phosphosites_tsv,
                                                 pssm_dict_h5 )
        tsvgz = pssm_bg_scores.tsvgz
        h5 = pssm_bg_scores.h5

    emit:
        tsvgz
        h5

}


workflow Y_PSSM_BACKGROUND_SCORES {

    take:
        human_phosphosites_tsv
        pssm_dict_h5

    main:
        pssm_bg_scores = y_pssm_background_scores( human_phosphosites_tsv,
                                                 pssm_dict_h5 )
        tsvgz = pssm_bg_scores.tsvgz
        h5 = pssm_bg_scores.h5

    emit:
        tsvgz
        h5

}

workflow PUBLISH_CONFIG {

    main:
        config_ch = Channel.fromPath("${projectDir}/nextflow.config")
        val_ch = Channel.of('workflow/nextflow.config')
        
        publish( config_ch, val_ch )

}

workflow {

    // get alias table from STRING
    string_id_dict = GET_STRING_ID_DICT()

    // get serine/threonine PSSMs
    ser_thr_kinases_pssm_dict_h5 = SER_THR_KINASES_PSSM( string_id_dict )

    // get tyrosine PSSMs
    tyr_kinases_pssm_dict_h5 = TYR_KINASES_PSSM( string_id_dict )

    // get human kinase-phosphosite associations from phosphositeplus
    psp_data = PSP( psp_kinase_substrate_gz,
                    psp_phosphosites_dataset_gz,
                    string_id_dict )
    psp_kin_sub_clusters = psp_data.kin_sub_clusters
    psp_human_phosphosites = psp_data.human_phosphosites

    // compute a priori pssm score distributions for Ser/Thr kinases
    s_t_pssm_bg_scores = S_T_PSSM_BACKGROUND_SCORES( psp_human_phosphosites,
                                                     ser_thr_kinases_pssm_dict_h5 )

    // compute a priori pssm score distributions for Tyr kinases
    y_pssm_bg_scores = Y_PSSM_BACKGROUND_SCORES( psp_human_phosphosites,
                                                 tyr_kinases_pssm_dict_h5 )

    PUBLISH_CONFIG()

}