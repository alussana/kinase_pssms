#!/usr/bin/env python3

import pandas as pd
import numpy as np
import h5py
import sys


"""
dict_tsv = "input/gene_synomym2gene_name_dict.tsv"
kinase_annotations_tsv = 'input/kinase_annotations.tsv'
s_t_pssm_h5 = 'input/S_T_PSSMs.h5'
y_pssm_h5 = 'input/Y_PSSMs.h5'
metadata_h5 = 'input/kinase_metadata.h5'
out_h5 = 'kinase_metadata/kinase_metadata_annotated.h5'
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


def hdf5_to_dict(hdf5_file_path):
    """
    Convert an HDF5 file's contents into a nested dictionary recursively,
    with special handling for string datasets.

    Args:
        hdf5_file_path (str): Path to the HDF5 file

    Returns:
        dict: Nested dictionary containing the HDF5 file's structure and data
    """
    def _convert_group(h5_group):
        result = {}
        for key, item in h5_group.items():
            # If item is a group, recursively convert it
            if isinstance(item, h5py.Group):
                result[key] = _convert_group(item)
            # If item is a dataset
            elif isinstance(item, h5py.Dataset):
                # Check if the dataset contains string type
                if h5py.check_string_dtype(item.dtype):
                    # Handle scalar string
                    if item.shape == ():
                        value = item[()]
                        result[key] = (
                            value.decode("utf-8")
                            if isinstance(value, bytes)
                            else str(value)
                        )
                    # Handle array of strings
                    else:
                        value = item[:]
                        result[key] = [
                            v.decode("utf-8") if isinstance(v, bytes) else str(v)
                            for v in value
                        ]
                else:
                    # Handle non-string scalar datasets
                    if item.shape == ():
                        result[key] = item[()]
                    # Handle non-string array datasets
                    else:
                        result[key] = item[:]
                    # Convert numpy types to Python native types if possible
                    if hasattr(result[key], "tolist"):
                        result[key] = result[key].tolist()
        return result
    # Open the HDF5 file and convert its contents
    try:
        with h5py.File(hdf5_file_path, "r") as f:
            return _convert_group(f)
    except Exception as e:
        raise Exception(f"Error reading HDF5 file: {str(e)}")


"""def map_kinase_annotations_to_pssms(pssm_h5, annot_df, dict_df, default_specificity=np.nan, default_family=np.nan):
    pssms = h5py.File(pssm_h5, "r")
    kinase_names = list(pssms.keys())
    annot_kinase_names = list(annot_df["Gene name"])
    new_dict = {
        "specificity": {},
        "family": {},
    }
    for kinase_name in kinase_names:
        meta_kinase_name = None
        if kinase_name in annot_kinase_names:
            meta_kinase_name = kinase_name
        else:
            tr_list = list(dict_df.loc[dict_df["Gene name"]==kinase_name, "Gene synonym"].values)
            tr_list = tr_list + list(dict_df.loc[dict_df["Gene synonym"]==kinase_name, "Gene name"].values)
            for tr_str in tr_list:
                if tr_str in annot_kinase_names:
                    meta_kinase_name = tr_str
                    break
        if meta_kinase_name is None:
            new_dict["specificity"][kinase_name] = default_specificity
            new_dict["family"][kinase_name] = default_family
        else:
            var = annot_df.loc[annot_df["Gene name"]==meta_kinase_name, "Specificity"].values
            if len(var) > 0:
                new_dict["specificity"][kinase_name] = var[0]
            else:
                new_dict["specificity"][kinase_name] = default_specificity
            var = annot_df.loc[annot_df["Gene name"]==meta_kinase_name, "Family"].values
            if len(var) > 0:
                new_dict["family"][kinase_name] = var[0]
            else:
                new_dict["family"][kinase_name] = default_family
    return new_dict"""


def map_kinase_annotations_to_pssms(s_t_pssm_h5, y_pssm_h5, annot_df, dict_df, default_specificity=np.nan, default_family=np.nan):
    s_t_pssms = h5py.File(s_t_pssm_h5, "r")
    y_pssms = h5py.File(y_pssm_h5, "r")
    pssm_kinase_names = set(s_t_pssms.keys()).union((y_pssms.keys()))
    annot_kinase_names = set(annot_df["Gene name"])
    new_dict = {
        "specificity": {},
        "family": {},
    }
    for annot_kinase_name in annot_kinase_names:
        pssm_kinase_name = None
        if annot_kinase_name in pssm_kinase_names:
            pssm_kinase_name = annot_kinase_name
            pssm_kinase_names.remove(pssm_kinase_name)
        else:
            tr_list = list(dict_df.loc[dict_df["Gene name"]==annot_kinase_name, "Gene synonym"].values)
            tr_list = tr_list + list(dict_df.loc[dict_df["Gene synonym"]==annot_kinase_name, "Gene name"].values)
            for tr_str in tr_list:
                if tr_str in pssm_kinase_names:
                    pssm_kinase_name = tr_str
                    pssm_kinase_names.remove(pssm_kinase_name)
                    break
        if pssm_kinase_name is None:
            new_dict["specificity"][annot_kinase_name] = annot_df.loc[annot_df["Gene name"]==annot_kinase_name, "Specificity"].values[0]
            new_dict["family"][annot_kinase_name] = annot_df.loc[annot_df["Gene name"]==annot_kinase_name, "Family"].values[0]
        else:
            var = annot_df.loc[annot_df["Gene name"]==annot_kinase_name, "Specificity"].values
            if len(var) > 0:
                new_dict["specificity"][pssm_kinase_name] = var[0]
            else:
                new_dict["specificity"][pssm_kinase_name] = default_specificity
            var = annot_df.loc[annot_df["Gene name"]==annot_kinase_name, "Family"].values
            if len(var) > 0:
                new_dict["family"][pssm_kinase_name] = var[0]
            else:
                new_dict["family"][pssm_kinase_name] = default_family
    return new_dict



def main():
    dict_tsv = sys.argv[1]
    kinase_annotations_tsv = sys.argv[2]
    s_t_pssm_h5 = sys.argv[3]
    y_pssm_h5 = sys.argv[4]
    metadata_h5 = sys.argv[5]
    out_h5 = sys.argv[6]

    dict_df = pd.read_csv(dict_tsv, sep="\t", header=None)
    dict_df.columns = ["Gene synonym", "UniProt AC", "Gene name"]

    annot_df = pd.read_csv(kinase_annotations_tsv, sep="\t", header=None)
    annot_df.columns = ["Gene name", "Gene name tr from AC", "Specificity", "Family"]

    kinase_metadata_dict = hdf5_to_dict(metadata_h5)

    annotations = map_kinase_annotations_to_pssms(s_t_pssm_h5, y_pssm_h5, annot_df, dict_df)
    
    kinase_metadata_dict["specificity"] = {**annotations["specificity"]}
    kinase_metadata_dict["family"] = {**annotations["family"]}
    
    # Save the dictionary to an HDF5 file
    with h5py.File(out_h5, 'w') as h5file:
        save_dict_to_hdf5(h5file, '', kinase_metadata_dict)


if __name__ == "__main__":
    main()
