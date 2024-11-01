#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
get the collection of phosphosites in the human phosphoproteome from 
PhosphositePlus

filter for only human, non-ambigous phosphorylation sites

transform sequences to be all-uppercases

cut the information to only include the following:
.META:
1   UniProtAC_RESIDUE   P31946_T2
2   Sequence (-5,4)     _____MTMDK
*/
process hs_phosphoproteome {

    publishDir "${out_dir}",
            pattern: 'datasets/psp/*.tsv',
            mode: 'copy'

    input:
        path 'input/Phosphorylation_site_dataset.gz'

    output:
        path 'datasets/psp/hs_phosphoproteome.tsv'

    script:
    """
    mkdir -p datasets/psp

    zcat input/Phosphorylation_site_dataset.gz \
        | sed '1,4d' \
        | cut -f3,5,7,10,15 \
        | grep "human" \
        | awk -F '\\t' '\$5==0'  \
        | sed 's/-p//' \
        | cut -f-2,4 \
        | sed 's/\\t/_/' \
        | awk '{print \$1"\\t"toupper(substr(\$2,3,10))}' \
        > datasets/psp/hs_phosphoproteome.tsv
    """

}