#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
process experimentally validated human kinase-substrates associations from
PhosphositePlus

Filter for only human kinases and substrates

Remove KIN_ORGANISM and SUB_ORGANISM fields

.META:
 1  GENE
 2  KINASE
 3  KIN_ACC_ID
 4  SUBSTRATE
 5  SUB_GENE_ID
 6  SUB_ACC_ID
 7  SUB_GENE
 8  SUB_MOD_RSD
 9  SITE_GRP_ID
10  SITE_+/-7_AA
11  DOMAIN
12  IN_VIVO_RXN
13  IN_VITRO_RXN
14  CST_CAT#
*/
process psp_kinase_substrates {

    publishDir "${out_dir}",
            pattern: 'datasets/psp/*.tsv',
            mode: 'copy'

    input:
        path 'input/psp_kinase_substrates.gz'

    output:
        path 'datasets/psp/kinase_substrates.tsv'

    script:
    """
    mkdir -p datasets/psp

    cat \
        <(zcat input/psp_kinase_substrates.gz \
            | sed '1,3d' | sed -n '1p' \
            | cut -f-3,5-8,10- ) \
        <(zcat input/psp_kinase_substrates.gz \
            | awk '{if(\$4=="human" && \$9=="human"){print \$0}}' \
            | cut -f-3,5-8,10- ) \
        > datasets/psp/kinase_substrates.tsv
    """

}

/*
.META:
1   Kinase Name
2   tab-separated list of substrates (UniProtAC_ResiduePosition)
*/
process kin_phos_clusters {

    publishDir "${out_dir}",
            pattern: 'datasets/psp/*.tsv',
            mode: 'copy'

    input:
        path 'input/kinase_substrates.tsv'

    output:
        path 'datasets/psp/kin_phos_clusters.tsv'

    script:
    """
    mkdir -p datasets/psp/

    cat input/kinase_substrates.tsv \
        | sed '1d' \
        | cut -f2 \
        | sort | uniq \
        > kinase_names.txt
    
    cat input/kinase_substrates.tsv \
        | sed '1d' \
        | cut -f2,6,8 \
        | awk '{print \$1"\t"\$2"_"\$3}' \
        > kinase_substrates_tab.tsv

    for kinase in \$(cat kinase_names.txt); do \
        sub_list=\$(cat kinase_substrates_tab.tsv \
            | grep "\$kinase" \
            | cut -f2 \
            | tr '\\n' '\\t'); \
        echo -e "\${kinase}\\t\${sub_list}" >> datasets/psp/kin_phos_clusters.tsv; \
    done

    """

}

/*
Translate gene synonym into gene names for kinases

.META:
1   Kinase Name
2   tab-separated list of substrates (UniProtAC_ResiduePosition)
*/
process translate_kin_phos_clusters_w_uniprot {

    publishDir "${out_dir}",
            pattern: 'datasets/psp/*.tsv',
            mode: 'copy'

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'datasets/psp/kin_phos_clusters_tr.tsv'

    script:
    """
    mkdir -p datasets/psp/

    cat input/dict.tsv \
        | grep -v PAK1 \
        | grep -v PDHK1 \
        > dict.tsv

    translator.py \
        dict.tsv \
        input/file.tsv \
        1 \
        3 \
        1 \
        1 \
        > datasets/psp/kin_phos_clusters_tr.tsv

    """

}


/*
Translate kinase aliases with Ensembl HGNC gene symbol

.META:
1   Kinase Name
2   tab-separated list of substrates (UniProtAC_ResiduePosition)
*/
process translate_kin_phos_clusters {

    publishDir "${out_dir}",
            pattern: 'datasets/psp/*.tsv',
            mode: 'copy'

    input:
        path 'input/k_p_clusters.tsv'
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path 'datasets/psp/kin_phos_clusters_tr.tsv'

    script:
    """
    mkdir -p datasets/psp/

    translate_k_p_clusters.py \
        > datasets/psp/kin_phos_clusters_tr.tsv \
        2> kinases_w_exceptions.txt
    """

}


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