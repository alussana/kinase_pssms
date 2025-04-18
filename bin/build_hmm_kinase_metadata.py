#!/usr/bin/env python3

import pandas as pd
import sys


"""
kinases_tsv = 'input/kinases.tsv'
aloops_tsv = "kinase_domains_hmm_aloops.tsv"
"""


def main():
    kinases_tsv = sys.argv[1]
    aloops_tsv = sys.argv[2]

    df = pd.read_csv(kinases_tsv, sep="\t")
    aloops_df = pd.read_csv(aloops_tsv, sep="\t", header=None)
    aloops_df.columns = ["Entry", "A-loop HMM"]

    df = df.merge(aloops_df, on="Entry", how="outer")

    df["A-loop"] = df["A-loop HMM"].combine_first(df["A-loop"])
    
    df = df.drop(columns=["A-loop HMM"])

    print(df.to_csv(sep="\t", index=False, header=True))


if __name__ == "__main__":
    main()
