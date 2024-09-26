#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
.META:
1. UniProt_AC 
2. ID_type 
3. ID
*/
process dl_human {

    publishDir "${out_dir}",
                pattern: 'datasets/uniprot/HUMAN_9606_idmapping.dat.gz',
                mode: 'copy'

    output:
        path 'datasets/uniprot/HUMAN_9606_idmapping.dat.gz'

    shell:
    """
    mkdir -p datasets/uniprot

    wget -P datasets/uniprot ${params.url_gene_id_dict_human}
    """

}

/*
Create dictionary UniProc AC --> Ensembl Protein

ENSP are truncated to remove the version number

.META:
1. UniProt_AC
2. ENSP
*/
process uniprotac_to_ensp_dict {

    publishDir "${out_dir}",
                pattern: 'datasets/uniprot/UniProtAC2ENSP_dict.tsv',
                mode: 'copy'

    input:
        path 'input/HUMAN_9606_idmapping.dat.gz'

    output:
        path 'datasets/uniprot/UniProtAC2ENSP_dict.tsv'

    shell:
    """
    mkdir -p datasets/uniprot

    zcat input/HUMAN_9606_idmapping.dat.gz \
        | awk '\$2=="Ensembl_PRO"{print \$1"\\t"\$3}' \
        | sed 's/\\..*\$//' \
        > datasets/uniprot/UniProtAC2ENSP_dict.tsv
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
process IDa2uniprot2IDb {

    publishDir "${out_dir}",
                pattern: "uniprot/${IDa}2uniprot2${IDb}.tsv",
                mode: 'copy'

    input:
        path 'input/mapping.tsv.gz'
        val IDa
        val IDb

    output:
        path "datasets/uniprot/${IDa}2uniprot2${IDb}.tsv"

    script:
    """
    mkdir -p datasets/uniprot

    zcat input/mapping.tsv.gz \
        | grep -w -f <(echo -e "${IDa}\n${IDb}") \
        | gzip > mapping.tsv.gz

    IDa2uniprot2IDb.py \
        mapping.tsv.gz \
        ${IDa} \
        ${IDb} \
        > datasets/uniprot/${IDa}2uniprot2${IDb}.tsv
    """

}