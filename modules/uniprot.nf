#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
.META:
1. UniProtKB-AC 
2. ID_type 
3. ID
*/
process dl_human {

    publishDir "${out_dir}",
                pattern: 'uniprot/HUMAN_9606_idmapping.dat.gz',
                mode: 'copy'

    output:
        path 'uniprot/HUMAN_9606_idmapping.dat.gz'

    shell:
    """
    mkdir -p uniprot

    wget -P uniprot ${params.human_url}
    """

}


/*
Use the UniProt API to retrieve a list of protein AC for reviewed entries
belonging to the human proteome
*/
process dl_ref_proteome_ac_list {

    publishDir "${out_dir}",
            pattern: 'databases/uniprot/hs_proteome_ac_list.txt',
            mode: 'copy'

    output:
        path 'databases/uniprot/hs_proteome_ac_list.txt'

    shell:
    """
    mkdir -p databases/uniprot

    wget -O databases/uniprot/hs_proteome_ac_list.txt '${params.human_reviewed_proteome_query}'
    """

}


/*
3-fields tab-delimited file

provides mapping in this order:

IDa --->  UniProt AC ---> IDb

all ids of nomenclature IDa are mapped to the corresponding UniProt ACs,
and the UniProt ACs are then mapped to the corresponding ids of nomenclature
IDb. For example, all ENSP ids are mapped to their UniProt AC. Then, each of
the UniProt AC identifiers is mapped to HGNC ids.

.META:
1. IDa
2. UniProt AC referred to IDa
3. IDb referred to UniProt AC
*/
process IDa2uniprot_ref2IDb {

    publishDir "${out_dir}",
                pattern: "uniprot/${IDa}2uniprot_ref2${IDb}.tsv",
                mode: 'copy'

    input:
        path 'input/mapping.tsv.gz'
        val IDa
        val IDb
        path 'input/hs_proteome_ac_list.txt'

    output:
        path "uniprot/${IDa}2uniprot_ref2${IDb}.tsv"

    script:
    """
    mkdir -p uniprot

    zcat input/mapping.tsv.gz \
        | grep -w -f <(echo -e "${IDa}\n${IDb}") \
        | gzip > mapping.tsv.gz

    IDa2uniprot2IDb.py \
        mapping.tsv.gz \
        ${IDa} \
        ${IDb} \
        | grep -w -f input/hs_proteome_ac_list.txt \
        > uniprot/${IDa}2uniprot_ref2${IDb}.tsv
    """

}


/*
filter the UniProt id mapping information to only incude mapping to:

Gene_Name
Gene_Synonym
*/
process uniprot2gene_name_and_synonym {

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/uniprot2gene_synonym.tsv',
                mode: 'copy'

    input:
        path 'input/HUMAN_9606_idmapping.dat.gz'

    output:
        path 'databases/uniprot/uniprot2gene_synonym.tsv'

    script:
    """
    mkdir -p databases/uniprot

    echo -e "Gene_Name\\nGene_Synonym" \
        > idtypes.txt

    zcat input/HUMAN_9606_idmapping.dat.gz \
        | grep -w -f idtypes.txt \
        | sed -r 's/(UniProtKB-ID\\t[0-9A-Z]+)_HUMAN\$/\\1/' \
        > databases/uniprot/uniprot2gene_synonym.tsv
    """

}


/*
translate all the words in the first field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the third column

the second column of input/dict.tsv specifies whether the translation is to a
"Gene_Name" or a "Gene_Synonym"

UniProt AC ids are preferentially translated to Gene_Name; if a gene name is
not found in the dictionary, then we look for a Gene_Synonym 

keep untranslated rows
*/
process translate_ac_to_gene_name {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    translator_ac_2_gene_name_or_synonym.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        3 \
        1 \
        1 \
        2 \
        > translated_file.tsv
    """

}


/*
translate all the words in the second field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the third column

the second column of input/dict.tsv specifies whether the translation is to a
"Gene_Name" or a "Gene_Synonym"

UniProt AC ids are preferentially translated to Gene_Name; if a gene name is
not found in the dictionary, then we look for a Gene_Synonym 

keep untranslated rows
*/
process translate_ac_to_gene_name_2 {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    translator_ac_2_gene_name_or_synonym.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        3 \
        2 \
        1 \
        2 \
        > translated_file.tsv
    """

}