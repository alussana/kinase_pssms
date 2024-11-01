#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_ser_thr_kinome_2023_suppl_table_2 } from './modules/pssm'
include { get_ser_thr_kinase_pssm } from './modules/pssm'
include { get_tyr_kinome_2024_suppl_table_2 } from './modules/pssm'
include { get_tyr_kinase_pssm } from './modules/pssm'

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

    main:
        kinome_2023_suppl_table_2 = get_ser_thr_kinome_2023_suppl_table_2()
        pssms = get_ser_thr_kinase_pssm( kinome_2023_suppl_table_2 )
                    .pssm_dict_h5

    emit:
        pssms

}


workflow TYR_KINASES_PSSM {

    main:
        kinome_2024_suppl_table_2 = get_tyr_kinome_2024_suppl_table_2()
        pssms = get_tyr_kinase_pssm( kinome_2024_suppl_table_2 )
                    .pssm_dict_h5

    emit:
        pssms

}


workflow PSP {

    take:
        psp_phosphosites_dataset_gz

    main:
        human_phosphosites = hs_phosphoproteome( psp_phosphosites_dataset_gz )

    emit:
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

    // get serine/threonine PSSMs
    ser_thr_kinases_pssm_dict_h5 = SER_THR_KINASES_PSSM()


    // get tyrosine PSSMs
    tyr_kinases_pssm_dict_h5 = TYR_KINASES_PSSM()


    // get human phosphosites from phosphositeplus
    human_phosphosites = PSP( psp_phosphosites_dataset_gz )


    // compute a priori pssm score distributions for Ser/Thr kinases
    s_t_pssm_bg_scores = S_T_PSSM_BACKGROUND_SCORES( human_phosphosites,
                                                     ser_thr_kinases_pssm_dict_h5 )


    // compute a priori pssm score distributions for Tyr kinases
    y_pssm_bg_scores = Y_PSSM_BACKGROUND_SCORES( human_phosphosites,
                                                 tyr_kinases_pssm_dict_h5 )


    PUBLISH_CONFIG()

}