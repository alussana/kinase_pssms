#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
.META:
1   #string_protein_id  9606.ENSP00000000233
2   alias               2B6H
3   source              Ensembl_PDB

=== Supported sources ===
$ zcat datasets/string/9606.protein.aliases.v12.0.txt.gz \
    | sed '1d' | cut -f3 | sort | uniq \
    > string_alias_sources.txt
$ cat string_alias_sources.txt
BioMart_HUGO
BioMart_HUGO_Paralog
Ensembl
Ensembl_archive_transcript
Ensembl_archive_translation
Ensembl_ArrayExpress
Ensembl_CCDS
Ensembl_ChEMBL
Ensembl_Clone_based_ensembl_gene
Ensembl_Clone_based_ensembl_transcript
Ensembl_DBASS3
Ensembl_DBASS5
Ensembl_description
Ensembl_EMBL
Ensembl_ENS_LRG_gene
Ensembl_ENS_LRG_transcript
Ensembl_EntrezGene
Ensembl_EntrezGene_Paralog
Ensembl_EntrezGene_trans_name
Ensembl_external_synonym_EntrezGene
Ensembl_external_synonym_HGNC
Ensembl_gene
Ensembl_GeneDB
Ensembl_HGNC
Ensembl_HGNC_alias_name
Ensembl_HGNC_alias_symbol
Ensembl_HGNC_ccds_id
Ensembl_HGNC_ena
Ensembl_HGNC_ensembl_gene_id
Ensembl_HGNC_entrez_id
Ensembl_HGNC_enzyme_id
Ensembl_HGNC_hgnc_id
Ensembl_HGNC_merops
Ensembl_HGNC_name
Ensembl_HGNC_omim_id
Ensembl_HGNC_prev_name
Ensembl_HGNC_prev_symbol
Ensembl_HGNC_pseudogene.org
Ensembl_HGNC_refseq_accession
Ensembl_HGNC_symbol
Ensembl_HGNC_trans_name
Ensembl_HGNC_ucsc_id
Ensembl_HGNC_uniprot_ids
Ensembl_HGNC_vega_id
Ensembl_HPA
Ensembl_MEROPS
Ensembl_PDB
Ensembl_protein_id
Ensembl_Reactome
Ensembl_Reactome_gene
Ensembl_Reactome_transcript
Ensembl_RefSeq
Ensembl_RefSeq_short
Ensembl_RNAcentral
Ensembl_transcript
Ensembl_translation
Ensembl_UCSC
Ensembl_UniParc
Ensembl_UniProt
Ensembl_WikiGene
KEGG_EC
KEGG_GENEID
KEGG_KEGGID
KEGG_KEGGID_SHORT
KEGG_NAME
KEGG_NAME_Paralog
KEGG_NAME_SYNONYM
KEGG_NCBI
KEGG_PRODUCT
UniProt_AC
UniProt_DE_AltName_Allergen
UniProt_DE_AltName_CD_antigen
UniProt_DE_AltName_EC
UniProt_DE_AltName_Full
UniProt_DE_AltName_INN
UniProt_DE_AltName_Short
UniProt_DE_RecName_EC
UniProt_DE_RecName_Full
UniProt_DE_RecName_Short
UniProt_DE_SubName_Full
UniProt_DR_GeneID
UniProt_DR_neXtProt
UniProt_DR_PDB
UniProt_DR_RefSeq
UniProt_DR_UCSC
UniProt_GN_Name
UniProt_GN_Name_Paralog
UniProt_GN_ORFNames
UniProt_GN_Synonyms
UniProt_ID
*/
process dl_string_aliases {

    publishDir "${out_dir}",
                pattern: 'datasets/string/9606.protein.aliases.v12.0.txt.gz',
                mode: 'copy'

    output:
        path 'datasets/string/9606.protein.aliases.v12.0.txt.gz'

    shell:
    """
    mkdir -p datasets/string

    wget -P datasets/string ${params.url_string_id_map_human}
    """

}