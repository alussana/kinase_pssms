#!/usr/bin/env python3

import pandas as pd
import sys


"""
kinases_tsv = 'input/kinases.tsv'
"""


def main():
    kinases_tsv = sys.argv[1]

    df = pd.read_csv(kinases_tsv, sep="\t")

    entries = df["Entry"].to_list()
    sequences = df["Sequence"].to_list()

    for i in range(len(entries)):
        print(f">{entries[i]}")
        print(sequences[i])


if __name__ == "__main__":
    main()
