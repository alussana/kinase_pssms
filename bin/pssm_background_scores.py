#!/usr/bin/env python3

import pandas as pd
import numpy as np
import h5py
import sys
from multiprocessing import Pool

AA_LIST = [
    "G",
    "P",
    "A",
    "V",
    "L",
    "I",
    "M",
    "C",
    "F",
    "Y",
    "W",
    "H",
    "K",
    "R",
    "Q",
    "N",
    "E",
    "D",
    "S",
    "T",
    "s",
    "t",
    "y",
]
POSITIONS_LIST = list(range(-5, 5))


def score_sequence(seq_str: str, pssm_df: pd.DataFrame):
    if len(seq_str) != len(pssm_df):
        raise Exception("Sequence length cannot be different from pssm length")
    """
    # transform pssm into probabilities, score sequence by multiplying values and scale by prob of random peptide
    pssm_df = pssm_df.apply(lambda x: x / x.sum(), axis=1)
    n_aa = len(AA_LIST)
    n_pos = len(seq_str)
    n_res = 0
    p = 1
    for i in range(n_pos):
        if seq_str[i] != '_':
            n_res = n_res + 1
            pos = list(pssm_df.index)[i]
            p = p * pssm_df.loc[pos, seq_str[i]] * 20
    """
    # score sequence by multiplying values and scale by prob of random peptide as described in Supplementary Note 2 of https://doi.org/10.1038/s41586-022-05575-3
    n_pos = len(seq_str)
    p = 1
    for i in range(n_pos):
        if seq_str[i] != "_":
            pos = list(pssm_df.index)[i]
            p = p * pssm_df.loc[pos, seq_str[i]]
    """
    # add scores instead of multiplying
    p = 0
    if pssm_df.loc[0, seq_str[5]] != 0:
        for i in range(len(seq_str)):
            if seq_str[i] != '_':
                pos = list(pssm_df.index)[i]
                p = p + pssm_df.loc[pos, seq_str[i]]
    """
    return p


def pssm_scoring(seq: str, pssm_df_dict: dict):
    record = {}
    for kinase in pssm_df_dict.keys():
        p = score_sequence(seq, pssm_df_dict[kinase])
        record[kinase] = p
    out_series = pd.Series(record)
    return out_series


def pssm_10000_quantiles(x: pd.Series):
    x.loc[x == 0] = np.nan
    x = x.dropna()
    quantiles_list = []
    for i in range(1, 10001):
        q = i / 10000
        quantiles_list.append(x.quantile(q))
    out_df = pd.Series(quantiles_list)
    out_df.name = x.name
    return out_df


def main():
    phosphosites_tsv = sys.argv[1]
    pssms_h5_file = sys.argv[2]
    n_proc = int(sys.argv[3])
    out_h5 = sys.argv[4]
    """
    phosphosites_tsv = 'input/phosphosites.tsv'
    pssms_h5_file = 'input/pssm_dict.h5'
    n_proc = 2
    out_h5 = 'PhosSEA/pssm_background_scores.h5'
    """
    pssms_h5 = h5py.File(pssms_h5_file, "r")
    pssm_df_dict = {}
    for kinase in pssms_h5.keys():
        pssm_df_dict[kinase] = pd.DataFrame(pssms_h5[kinase])
        pssm_df_dict[kinase].columns = AA_LIST
        pssm_df_dict[kinase].index = POSITIONS_LIST
    phosphosites = pd.read_csv(phosphosites_tsv, sep="\t", header=None, index_col=0)
    phosphosites = phosphosites[1]

    # score phosphosite sequences with each PSSM
    arg1 = [phosphosites[i] for i in range(len(phosphosites))]
    arg2 = [pssm_df_dict for i in range(len(phosphosites))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(pssm_scoring, zip(arg1, arg2))
    pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    pssm_scoring_df.index = list(range(len(phosphosites)))

    # print(pssm_scoring_df.to_csv(sep='\t', index=False, header=True))

    # create quantiles dataframe
    arg1 = [pssm_scoring_df[kinase] for kinase in list(pssm_scoring_df.columns)]
    with Pool(processes=n_proc) as pool:
        pssm_scores_quantiles_list = pool.starmap(pssm_10000_quantiles, zip(arg1))
    pssm_scores_quantiles_df = pd.concat(pssm_scores_quantiles_list, axis=1)

    print(pssm_scores_quantiles_df.to_csv(sep="\t", index=False, header=True))

    # Save the DataFrame to HDF5
    pssm_scores_quantiles_df.to_hdf(
        out_h5,
        key="pssm_scores",
        mode="w",
        format="table",
        complib="blosc",
        complevel=9,
    )


if __name__ == "__main__":
    main()
