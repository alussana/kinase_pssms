#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
Download Supplementary Table 2 from

Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature 613, 759–766 (2023).

<https://doi.org/10.1038/s41586-022-05575-3>
*/
process get_ser_thr_kinome_2023_suppl_table_2 {

    publishDir "${out_dir}", pattern: "datasets/41586_2022_5575_MOESM4_ESM.xlsx", mode: 'copy'

    output:
        path 'datasets/41586_2022_5575_MOESM4_ESM.xlsx'

    script:
    """
    mkdir -p datasets
    
    wget -P datasets/ "${params.url_ser_thr_kinome_2023_suppl_table_2}"
    """

}


/*
Parse output from get_ser_thr_kinome_2023_suppl_table_2()
*/
process get_ser_thr_kinase_pssm {

    publishDir "${out_dir}", pattern: "datasets/PSSMsm/*tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "datasets/PSSMs/S_T_PSSMs.h5", mode: 'copy'
    publishDir "${out_dir}", pattern: "datasets/PSSMs/logos/*.pdf", mode: 'copy'

    input:
        path 'input/kinome_2023_suppl_table_2.xlsx'

    output:
        path 'datasets/PSSMs/S_T_PSSMs.h5', emit: pssm_dict_h5
        path 'datasets/PSSMs/*.tsv'
        path 'datasets/PSSMs/logos/*.pdf'

    script:
    """
    mkdir -p datasets/PSSMs/logos

    get_pssm_ser_thr_kinome_2023_suppl_table_2.py
    """

}


/*
Download Supplementary Table 2 from

Yaron-Barir, T.M., Joughin, B.A., Huntsman, E.M. et al. 
The intrinsic substrate specificity of the human tyrosine kinome.
Nature 629, 1174–1181 (2024).

<https://doi.org/10.1038/s41586-024-07407-y>
*/
process get_tyr_kinome_2024_suppl_table_2 {

    publishDir "${out_dir}", pattern: "datasets/suppl_table_2.xlsx", mode: 'copy'

    output:
        path 'datasets/suppl_table_2.xlsx'

    script:
    """
    mkdir -p datasets
    
    wget -O datasets/suppl_table_2.xlsx "${params.url_tyr_kinome_2024_suppl_table_2}"
    """

}

/*
Parse PSSMs from the outpur of get_tyr_kinome_2024_suppl_table_2()
*/
process get_tyr_kinase_pssm {

    publishDir "${out_dir}", pattern: "datasets/PSSMs/Y_PSSMs.h5", mode: 'copy'
    publishDir "${out_dir}", pattern: "datasets/PSSMs/*tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "datasets/PSSMs/logos/*.pdf", mode: 'copy'

    input:
        path 'input/kinome_2024_suppl_table_2.xlsx'

    output:
        path 'datasets/PSSMs/Y_PSSMs.h5', emit: pssm_dict_h5
        path 'datasets/PSSMs/*.tsv'
        path 'datasets/PSSMs/logos/*.pdf'

    script:
    """
    CACHEBUST=1

    mkdir -p datasets/PSSMs/logos

    cp \$(readlink input/kinome_2024_suppl_table_2.xlsx) kinome_2024_suppl_table_2.xlsx

    get_pssm_tyr_kinome_2024_suppl_table_2.py
    """

}


/*
Translate kinase aliases with Ensembl HGNC gene symbol
*/
process s_t_pssm_h5_alias_mapping {

    publishDir "${out_dir}", pattern: "pssm/*.h5"
    publishDir "${out_dir}", pattern: "pssm/logos/*.pdf"

    input:
        path 'input/S_T_PSSMs.h5'
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path 'pssm/S_T_PSSMs.h5', emit: pssm_dict_h5
        path 'pssm/logos/*pdf'

    script:
    """
    mkdir -p pssm/logos

    cp \$(readlink input/S_T_PSSMs.h5) pssm/S_T_PSSMs.h5

    pssm_h5_alias_mapping.py \
         pssm/S_T_PSSMs.h5 \
        > kinases_w_exceptions.txt
    """

}


/*
Translate kinase aliases with Ensembl HGNC gene symbol
*/
process y_pssm_h5_alias_mapping {

    publishDir "${out_dir}", pattern: "pssm/*.h5"
    publishDir "${out_dir}", pattern: "pssm/logos/*.pdf"

    input:
        path 'input/Y_PSSMs.h5'
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path 'pssm/Y_PSSMs.h5', emit: pssm_dict_h5
        path 'pssm/logos/*pdf'

    script:
    """
    mkdir -p pssm/logos

    cp \$(readlink input/Y_PSSMs.h5) pssm/Y_PSSMs.h5

    pssm_h5_alias_mapping.py \
         pssm/Y_PSSMs.h5 \
        > kinases_w_exceptions.txt
    """

}