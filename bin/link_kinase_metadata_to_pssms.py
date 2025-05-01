#!/usr/bin/env python3

import pandas as pd
import numpy as np
import h5py
import sys


"""
dict_tsv = "input/gene_synomym2gene_name_dict.tsv"
kinase_metadata_tsv = 'input/kinase_metadata.tsv'
s_t_pssm_h5 = 'input/S_T_PSSMs.h5'
y_pssm_h5 = 'input/Y_PSSMs.h5'
out_h5 = 'kinase_metadata/kinases.h5'
"""


def save_dict_to_hdf5(h5file, path, dic):
    for key, value in dic.items():
        current_path = f"{path}/{key}"
        if isinstance(value, dict):
            # If the value is a dictionary, create a group and recurse
            group = h5file.create_group(current_path)
            save_dict_to_hdf5(h5file, current_path, value)
        else:
            # If the value is not a dictionary, save it as a dataset
            h5file.create_dataset(current_path, data=value)


def map_kinase_metadata_to_pssms(pssm_h5, meta_df, dict_df):
    pssms = h5py.File(pssm_h5, "r")
    kinase_names = list(pssms.keys())
    new_dict = {
        "full_seq": {},
        "kinase_domain_seq": {},
        "aloop_seq": {},
    }
    for kinase_name in kinase_names:
        meta_kinase_name = None
        if kinase_name in list(meta_df["Entry"]):
            meta_kinase_name = kinase_name
        else:
            tr_list = list(dict_df.loc[dict_df["Gene name"]==kinase_name, "Gene synonym"].values)
            tr_list = tr_list + list(dict_df.loc[dict_df["Gene synonym"]==kinase_name, "Gene name"].values)
            for tr_str in tr_list:
                if tr_str in kinase_names:
                    meta_kinase_name = tr_str
                    break
        if meta_kinase_name is None:
            new_dict["full_seq"][kinase_name] = np.nan
            new_dict["kinase_domain_seq"][kinase_name] = np.nan
            new_dict["aloop_seq"][kinase_name] = np.nan
        else:
            var = meta_df.loc[meta_df["Entry"]==meta_kinase_name, "Sequence"].values
            if len(var) > 0:
                new_dict["full_seq"][kinase_name] = var[0]
            else:
                new_dict["full_seq"][kinase_name] = np.nan
            var = meta_df.loc[meta_df["Entry"]==meta_kinase_name, "Kinase domain"].values
            if len(var) > 0:
                new_dict["kinase_domain_seq"][kinase_name] = var[0]
            else:
                new_dict["kinase_domain_seq"][kinase_name] = np.nan
            var = meta_df.loc[meta_df["Entry"]==meta_kinase_name, "A-loop"].values
            if len(var) > 0:
                new_dict["aloop_seq"][kinase_name] = var[0]
            else:
                new_dict["aloop_seq"][kinase_name] = np.nan
    return new_dict


def main():
    dict_tsv = sys.argv[1]
    kinase_metadata_tsv = sys.argv[2]
    s_t_pssm_h5 = sys.argv[3]
    y_pssm_h5 = sys.argv[4]
    out_h5 = sys.argv[5]

    dict_df = pd.read_csv(dict_tsv, sep="\t", header=None)
    dict_df.columns = ["Gene synonym", "UniProt AC", "Gene name"]

    meta_df = pd.read_csv(kinase_metadata_tsv, sep="\t")

    s_t_metadata = map_kinase_metadata_to_pssms(s_t_pssm_h5, meta_df, dict_df)
    y_metadata = map_kinase_metadata_to_pssms(y_pssm_h5, meta_df, dict_df)
    
    metadata = {}
    metadata["full_seq"] = {**s_t_metadata["full_seq"], **y_metadata["full_seq"]}
    metadata["kinase_domain_seq"] = {**s_t_metadata["kinase_domain_seq"], **y_metadata["kinase_domain_seq"]}
    metadata["aloop_seq"] = {**s_t_metadata["aloop_seq"], **y_metadata["aloop_seq"]}
    

    # Save the dictionary to an HDF5 file
    with h5py.File(out_h5, 'w') as h5file:
        save_dict_to_hdf5(h5file, '', metadata)


if __name__ == "__main__":
    main()
