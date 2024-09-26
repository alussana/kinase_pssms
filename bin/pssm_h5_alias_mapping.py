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


def plot_sequence_logo(
    pssm_df: pd.DataFrame,
    file_path: str = "logo.svg",
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
    alias_map_txt_gz = "input/9606.protein.aliases.v12.0.txt.gz"
    pssms_h5_file = sys.argv[1]
    
    alias_map_df = pd.read_csv(alias_map_txt_gz, sep="\t")
    
    with h5py.File(pssms_h5_file, "r+") as hdf_file:
        for key in hdf_file.keys():
            key_alias_df = alias_map_df.loc[alias_map_df["alias"]==key]
            string_id_array = key_alias_df["#string_protein_id"].unique()
            for string_id_str in string_id_array:
                try:
                    map_df = alias_map_df.loc[alias_map_df["#string_protein_id"]==string_id_str]
                    map_df.loc[map_df["source"]=="Ensembl_HGNC_symbol",]
                    new_key = map_df.loc[map_df["source"]=="Ensembl_HGNC_symbol","alias"].values[0]
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
                    break
                except:
                    print(key)
                    next


if __name__ == "__main__":
    main()
