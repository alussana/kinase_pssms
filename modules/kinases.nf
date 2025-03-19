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
