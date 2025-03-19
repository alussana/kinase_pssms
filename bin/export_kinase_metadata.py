#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import h5py
import sys


"""
kinases_tsv_gz = 'input/kinases.tsv.gz'
out_h5 = 'kinase_metadata/kinases.h5'
"""

def main():
    kinases_tsv_gz = sys.argv[1]
    out_h5 = sys.argv[2]

    df = pd.read_csv(kinases_tsv_gz, sep="\t")

    # Save the DataFrame to HDF5
    #kinases_h5.to_hdf(
    #    out_h5,
    #    key="UniProtAC",
    #    mode="w",
    #    format="table",
    #    complib="blosc",
    #    complevel=9,
    #)


if __name__ == "__main__":
    main()
