#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import h5py
import logomaker
import matplotlib.pyplot as plt


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


POS_LIST = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4]


alias_map_txt_gz = "input/9606.protein.aliases.v12.0.txt.gz"

pssms_h5_file = sys.argv[1]

alias_map_df = pd.read_csv(alias_map_txt_gz, sep="\t")

Ensembl_HGNC_symbol_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_HGNC_symbol", "alias"
].unique()

Ensembl_HGNC_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_HGNC", "alias"
].unique()

Ensembl_UniProt_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_UniProt", "alias"
].unique()


def find_translation(
    key,
):
    # The following mappings would collide if not for the ad-hoc flow control that follows:
    # Ser/Thr:
    # HGK     MAP4K4
    # NIK     MAP4K4
    # PDHK1   PDK1
    # PDK1    PDK1
    # Tyr:
    # EPHA3   EPHA3
    # ETK     EPHA3
    if key == "PDHK1":
        return key
    elif key == "HGK":
        return key
    elif key == "ETK":
        return key
    key_parts = key.split("_")
    key = key_parts[0]
    if (
        key in Ensembl_HGNC_symbol_list
        or key in Ensembl_HGNC_list
        or key in Ensembl_UniProt_list
    ):
        if len(key_parts) > 1:
            new_key = f"{key}_TYR"
        else:
            new_key = key
        return new_key
    key_alias_df = alias_map_df.loc[alias_map_df["alias"] == key]
    string_id_array = key_alias_df["#string_protein_id"].unique()
    if len(string_id_array) == 0:
        return key
    else:
        for string_id_str in string_id_array:
            try:
                map_df = alias_map_df.loc[
                    alias_map_df["#string_protein_id"] == string_id_str
                ]
                new_key = map_df.loc[map_df["source"] == "Ensembl_HGNC_symbol", "alias"]
                if len(new_key) == 0:
                    new_key = map_df.loc[map_df["source"] == "Ensembl_HGNC", "alias"]
                if len(new_key) == 0:
                    new_key = map_df.loc[map_df["source"] == "Ensembl_UniProt", "alias"]
                if len(new_key) == 0:
                    next
                new_key = new_key.values[0]
                if len(key_parts) > 1:
                    new_key = f"{new_key}_TYR"
                break
            except Exception as e:
                print(f"{key}: {e}")
                next
        return new_key


def plot_sequence_logo(
    pssm_df: pd.DataFrame,
    file_path: str = "logo.pdf",
    ylabel_str: str = "",
    title_str: str = "",
) -> None:

    central_residue_df = pd.DataFrame([pssm_df.loc[0,]])
    pssm_df = pssm_df.drop(0, axis=0)
    pssm_median_shifted_df = pssm_df.apply(lambda x: np.log2(x / x.median()), axis=1)
    max_pos_height_float = (
        pssm_median_shifted_df[pssm_median_shifted_df > 0].T.sum().max()
    )
    central_residue_df = (
        central_residue_df * max_pos_height_float / central_residue_df.T.sum().values[0]
    )
    pssm_median_shifted_df = pd.concat(
        [pssm_median_shifted_df, central_residue_df], axis=0
    )
    pssm_median_shifted_df = pssm_median_shifted_df.sort_index()
    plt.clf()
    # create Logo object
    pssm_logo = logomaker.Logo(
        pssm_median_shifted_df,
        shade_below=0,
        fade_below=0,
        width=0.9,
        color_scheme="chemistry",
        figsize=[4, 2],
        flip_below=False,
        stack_order="big_on_top",
        center_values=False,
        show_spines=True,
        # font_name='Arial Rounded MT Bold')
    )
    # style using Logo methods
    pssm_logo.style_spines(visible=False)
    pssm_logo.style_spines(spines=["left", "bottom"], visible=True)
    pssm_logo.style_xticks(rotation=0, fmt="%d", anchor=0)
    # style using Axes methods
    pssm_logo.ax.set_ylabel(ylabel_str, labelpad=-1)
    pssm_logo.ax.set_title(title_str)
    pssm_logo.ax.xaxis.set_ticks_position("none")
    pssm_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig(file_path)
    plt.close()


def main():
    with h5py.File(pssms_h5_file, "r+") as hdf_file:
        for key in hdf_file.keys():
            new_key = find_translation(key)
            if new_key not in hdf_file:
                hdf_file.copy(key, new_key)
                del hdf_file[key]
            pssm_df = pd.DataFrame(hdf_file[new_key])
            pssm_df.columns = AA_LIST
            pssm_df.index = POS_LIST
            plot_sequence_logo(
                pssm_df,
                file_path=f"pssm/logos/{new_key}.pdf",
                ylabel_str="log2(Ratio to Median)",
                title_str=new_key,
            )
            print(f"{key}\t{new_key}")


if __name__ == "__main__":
    main()
