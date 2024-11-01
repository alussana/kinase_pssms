#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
compute kinase pssm scores for the human phosphoproteome
in order to obtain a background distribution of pssm scores for
each kinase; this data will be used in the kinase activity estimation
algorithm

pssm score of more than 200k human phosphosites are computed

for each kinase all non-zero scores are then used to extract the values of
10000 quantiles

these quantiles are used to represent the pssm score distributions 
*/
process s_t_pssm_background_scores {

    cpus "${params.n_cores}"

    memory "16G"

    publishDir "${out_dir}", pattern: "pssm_bg_scores/*.tsv.gz", mode: 'copy'
    publishDir "${out_dir}", pattern: "pssm_bg_scores/*.h5", mode: 'copy'

    input:
        path 'input/phosphosites.tsv'
        path 'input/pssm_dict.h5'

    output:
        path "pssm_bg_scores/*.tsv.gz", emit: tsvgz
        path "pssm_bg_scores/*.h5", emit: h5

    script:
    """
    mkdir -p pssm_bg_scores

    pssm_background_scores.py \
        input/phosphosites.tsv \
        input/pssm_dict.h5 \
        ${params.n_cores} \
        pssm_bg_scores/S_T_PSSMs_score_quantiles.h5 \
        | gzip > pssm_bg_scores/S_T_PSSM_score_quantiles.tsv.gz
    """

}


/*
compute kinase pssm scores for the human phosphoproteome
in order to obtain a background distribution of pssm scores for
each kinase; this data will be used in the kinase activity estimation
algorithm

pssm score of more than 200k human phosphosites are computed

for each kinase all non-zero scores are then used to extract the values of
10000 quantiles

these quantiles are used to represent the pssm score distributions 
*/
process y_pssm_background_scores {

    cpus "${params.n_cores}"

    memory "16G"

    publishDir "${out_dir}", pattern: "pssm_bg_scores/*.tsv.gz", mode: 'copy'
    publishDir "${out_dir}", pattern: "pssm_bg_scores/*.h5", mode: 'copy'

    input:
        path 'input/phosphosites.tsv'
        path 'input/pssm_dict.h5'

    output:
        path "pssm_bg_scores/*.tsv.gz", emit: tsvgz
        path "pssm_bg_scores/*.h5", emit: h5

    script:
    """
    mkdir -p pssm_bg_scores

    pssm_background_scores.py \
        input/phosphosites.tsv \
        input/pssm_dict.h5 \
        ${params.n_cores} \
        pssm_bg_scores/Y_PSSMs_score_quantiles.h5 \
        | gzip > pssm_bg_scores/Y_PSSM_score_quantiles.tsv.gz
    """

}