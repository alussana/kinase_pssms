#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
Fetch human kinases information from UniProt

.META
1   Entry       O00141       
2   Gene Names  SGK1 SGK
3   Sequence    MTVKTEAAKGTLTYSRMRGMVAILIAFMKQRRMGLNDFIQKIANNSYACKHPEVQSILKISQPQEPELMNANPSPPPSPSQQINLGPSSNPHAKPSDFHFLKVIGKGSFGKVLLARHKAEEVFYAVKVLQKKAILKKKEEKHIMSERNVLLKNVKHPFLVGLHFSFQTADKLYFVLDYINGGELFYHLQRERCFLEPRARFYAAEIASALGYLHSLNIVYRDLKPENILLDSQGHIVLTDFGLCKENIEHNSTTSTFCGTPEYLAPEVLHKQPYDRTVDWWCLGAVLYEMLYGLPPFYSRNTAEMYDNILNKPLQLKPNITNSARHLLEGLLQKDRTKRLGAKDDFMEIKSHVFFSLINWDDLINKKITPPFNPNVSGPNDLRHFDPEFTEEPVPNSIGKSPDSVLVTASVKEAAEAFLGFSYAPPTDSFL
4   Domain[FT]  DOMAIN 98..355; /note="Protein kinase"; /evidence="ECO:0000255|PROSITE-ProRule:PRU00159"; DOMAIN 356..431; /note="AGC-kinase C-terminal"; /evidence="ECO:0000255|PROSITE-ProRule:PRU00618"
*/
process get_human_kinases {

    publishDir "${out_dir}", pattern: "datasets/uniprot/kinases.tsv.gz", mode: 'copy'

    output:
        path "datasets/uniprot/kinases.tsv.gz"

    script:
    """
    CACHEBUST=1

    mkdir -p datasets/uniprot
    
    wget -O datasets/uniprot/kinases.tsv.gz "${params.url_uniprot_kinases}"
    """

}


/*
Build kinase metadata

*/
process build_kinase_metadata {

    publishDir "${out_dir}", pattern: "kinase_metadata/*.tsv", mode: 'copy'

    input:
        path "input/kinases.tsv.gz"

    output:
        path "kinase_metadata/kinases.tsv"

    script:
    """
    mkdir -p kinase_metadata
    
    build_kinase_metadata.py input/kinases.tsv.gz > kinase_metadata/kinases.tsv
    """

}


/*
Build kinase metadata h5

ID-map the kinase metadata to the PSSMs
*/
process link_kinase_metadata_to_pssms {

    publishDir "${out_dir}", pattern: "kinase_metadata/*.h5", mode: 'copy'

    input:
        path "input/kinase_metadata.tsv"
        path "input/gene_synomym2gene_name_dict.tsv"
        path "input/S_T_PSSMs.h5"
        path "input/Y_PSSMs.h5"

    output:
        path "kinase_metadata/kinase_metadata.h5"

    script:
    """
    mkdir -p kinase_metadata
    
    link_kinase_metadata_to_pssms.py \
        input/gene_synomym2gene_name_dict.tsv \
        input/kinase_metadata.tsv \
        input/S_T_PSSMs.h5 \
        input/Y_PSSMs.h5 \
        kinase_metadata/kinase_metadata.h5
    """

}
