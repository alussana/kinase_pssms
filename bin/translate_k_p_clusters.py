#!/usr/bin/env python3

import pandas as pd
import sys


def main():
    alias_map_txt_gz = "input/9606.protein.aliases.v12.0.txt.gz"
    k_p_file = "input/k_p_clusters.tsv"

    alias_map_df = pd.read_csv(alias_map_txt_gz, sep="\t")

    with open(k_p_file, "r") as k_p_fh:
        for line in k_p_fh:
            fields_list = line.removesuffix("\n").split("\t")
            key = fields_list.pop(0)
            fields_str = "\t".join(fields_list)
            key_alias_df = alias_map_df.loc[alias_map_df["alias"] == key]
            string_id_array = key_alias_df["#string_protein_id"].unique()
            for string_id_str in string_id_array:
                try:
                    map_df = alias_map_df.loc[
                        alias_map_df["#string_protein_id"] == string_id_str
                    ]
                    map_df.loc[map_df["source"] == "Ensembl_HGNC_symbol",]
                    new_key = map_df.loc[
                        map_df["source"] == "Ensembl_HGNC_symbol", "alias"
                    ].values[0]
                    print(f"{new_key}\t{fields_str}")
                    break
                except:
                    print(f"{key}\t{fields_str}")
                    print(key, file=sys.stderr)
                    next


if __name__ == "__main__":
    main()
