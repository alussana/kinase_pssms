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
Fetch the PFAM kinase domain HMM (PF00069)
*/
process get_kinase_domain_hmm {

    publishDir "${out_dir}", pattern: "datasets/interpro/PF00069.hmm", mode: 'copy'

    output:
        path "datasets/interpro/PF00069.hmm"

    script:
    """
    mkdir -p datasets/interpro
    
    wget -O datasets/interpro/PF00069.hmm.gz "${params.url_pfam_prot_kinase_dom_hmm}"

    zcat datasets/interpro/PF00069.hmm.gz > datasets/interpro/PF00069.hmm
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


/*
[...]
*/
process make_kinase_sequences_fasta {

    publishDir "${out_dir}", pattern: "kinase_metadata/kinases.fa", mode: 'copy'

    input:
        path "input/kinases.tsv"

    output:
        path "kinase_metadata/kinases.fa"

    script:
    """
    mkdir -p kinase_metadata
    
    make_kinase_domains_fasta.py input/kinases.tsv > kinase_metadata/kinases.fa
    """

}


/*
Use HMMER to generate a multiple sequence alignment based on an HMM profile
*/
process run_hmmsearch {

    publishDir "${out_dir}", pattern: "kinase_metadata/kinase_domains.hmmsearch", mode: 'copy'

    input:
        path "input/kinase_domain.hmm"
        path "input/kinases.fa"

    output:
        path "kinase_metadata/kinase_domains.hmmsearch"

    script:
    """
    mkdir -p kinase_metadata
    
    hmmsearch --domtblout kinase_metadata/kinase_domains.hmmsearch input/kinase_domain.hmm input/kinases.fa > kinase_metadata/kinase_domains.hmmsearch
    """

}


/*
Retrieve kinase domain(s) HMM matches from original targets
using the domtblout hmmserach output
*/
process extract_kinase_domain_matches {

    publishDir "${out_dir}", pattern: "kinase_metadata/kinase_domains_hmm_matches.fa", mode: 'copy'

    input:
        path "input/kinase_domains.hmmsearch"
        path "input/kinases.fa"

    output:
        path "kinase_metadata/kinase_domains_hmm_matches.fa"

    script:
    """
    mkdir -p kinase_metadata
    
    extract_kinase_domain_matches.py \
        input/kinase_domains.hmmsearch \
        input/kinases.fa \
        kinase_metadata/kinase_domains_hmm_matches.fa
    """

}


/*
heuristically look for A-loops
*/
process extract_aloops {

    publishDir "${out_dir}", pattern: "kinase_metadata/kinase_domains_hmm_matches_aloops.fa", mode: 'copy'

    input:
        path "input/kinase_domains_hmm_matches.fa"

    output:
        path "kinase_metadata/kinase_domains_hmm_matches_aloops.fa"

    script:
    """
    mkdir -p kinase_metadata
    
    extract_aloops.py \
        input/kinase_domains_hmm_matches.fa \
        > kinase_metadata/kinase_domains_hmm_matches_aloops.fa
    """

}


/*
remove ambiguous A-loops, i.e. kinases with more than one A-loop are not
considered to avoid ambiguity

also remove instances where no A-loop has been found
*/
process remove_ambiguous_aloops {

    publishDir "${out_dir}", pattern: "kinase_metadata/kinase_domains_hmm_aloops.fa", mode: 'copy'

    input:
        path "input/kinase_domains_hmm_matches_aloops.fa"

    output:
        path "kinase_metadata/kinase_domains_hmm_aloops.fa"

    script:
    """
    mkdir -p kinase_metadata
    
    cat input/kinase_domains_hmm_matches_aloops.fa \
        | awk '{if (\$0 ~ /^>/) printf "%s ", \$0; else print \$0}' \
        | awk '\$2!="none"' \
        | cut -d "_" -f1 \
        | sort | uniq -d \
        > ambiguous_aloops_ac_list.txt
        
    cat input/kinase_domains_hmm_matches_aloops.fa \
        | grep -n -f ambiguous_aloops_ac_list.txt \
        | awk -v n=1 -F":" '{start = \$1; end = start + n; print start "," end "d"}' \
        | sed -E -f - input/kinase_domains_hmm_matches_aloops.fa \
        > kinase_metadata/kinase_domains_hmm_aloops_w_none.fa

    cat kinase_metadata/kinase_domains_hmm_aloops_w_none.fa \
        | grep -n "none" \
        | awk -v n=1 -F":" '{start = \$1 - n; end = \$1; print start "," end "d"}' \
        | sed -E -f - kinase_metadata/kinase_domains_hmm_aloops_w_none.fa \
        > kinase_metadata/kinase_domains_hmm_aloops.fa
    """

}


/*
Build kinase metadata
*/
process build_hmm_kinase_metadata {

    publishDir "${out_dir}", pattern: "kinase_metadata/*.tsv", mode: 'copy'

    input:
        path "input/kinases.tsv"
        path "input/kinase_domains_hmm_aloops.fa"

    output:
        path "kinase_metadata/hmm_kinases.tsv"

    script:
    """
    mkdir -p kinase_metadata
    
    cat input/kinase_domains_hmm_aloops.fa  \
        | awk '{if (\$0 ~ /^>/) printf "%s\\t", \$0; else print \$0}' \
        | tr -d '>' \
        | sed -e 's/_[0-9]*-[0-9]*//' \
        > kinase_domains_hmm_aloops.tsv
    
    build_hmm_kinase_metadata.py \
        input/kinases.tsv \
        kinase_domains_hmm_aloops.tsv \
        > kinase_metadata/hmm_kinases.tsv
    """

}


/*
Build kinase metadata h5

ID-map the kinase metadata to the PSSMs
*/
process link_hmm_kinase_metadata_to_pssms {

    publishDir "${out_dir}", pattern: "kinase_metadata/*.h5", mode: 'copy'

    input:
        path "input/kinase_metadata.tsv"
        path "input/gene_synomym2gene_name_dict.tsv"
        path "input/S_T_PSSMs.h5"
        path "input/Y_PSSMs.h5"

    output:
        path "kinase_metadata/hmm_kinase_metadata.h5"

    script:
    """
    mkdir -p kinase_metadata
    
    link_kinase_metadata_to_pssms.py \
        input/gene_synomym2gene_name_dict.tsv \
        input/kinase_metadata.tsv \
        input/S_T_PSSMs.h5 \
        input/Y_PSSMs.h5 \
        kinase_metadata/hmm_kinase_metadata.h5
    """

}


/*
Build kinase metadata h5 adding the annotations to an existing metadata h5 file

ID-map the kinase annotations to the PSSMs
*/
process link_kinase_annotations_to_pssms {

    publishDir "${out_dir}", pattern: "kinase_metadata/*.h5", mode: 'copy'

    input:
        path "input/kinase_annotations.tsv"
        path "input/gene_synomym2gene_name_dict.tsv"
        path "input/S_T_PSSMs.h5"
        path "input/Y_PSSMs.h5"
        path "input/kinase_metadata.h5"

    output:
        path "kinase_metadata/kinase_metadata_annotated.h5"

    script:
    """
    mkdir -p kinase_metadata
    
    link_kinase_annotations_to_pssms.py \
        input/gene_synomym2gene_name_dict.tsv \
        input/kinase_annotations.tsv \
        input/S_T_PSSMs.h5 \
        input/Y_PSSMs.h5 \
        input/kinase_metadata.h5 \
        kinase_metadata/kinase_metadata_annotated.h5
    """

}
