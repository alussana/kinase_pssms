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

def parse_domain_ft(s: str):
    if np.isnan(s):
        return None
    elif type(s) != str:
        raise TypeError(
            f"argument must be a string or nan, passed {type(s)}"
        )
    else:
        doms = s.split("DOMAIN ")
        doms.pop(0)
        features = []
        for d in doms:
            info = d.split("; ")
            info.pop()
            start, end = [int(n) for n in info[0].split("..")]
            name = info[1].split("\"")[1]
            features.append({
                "name": name,
                "start": start,
                "end": end,
            })
    return features
                    

def get_kinase_domain_seq(seq: str, domain_ft: str):
    if type(domain_ft) != str:
        return None
    else:
        doms = domain_ft.split("DOMAIN ")
        doms.pop(0)
        features = []
        for d in doms:
            info = d.split("; ")
            start, end = [int(n) for n in info[0].split("..")]
            name = info[1].split("\"")[1]
            features.append({
                "name": name,
                "start": start,
                "end": end,
            })
        # choose representative kinase domain
        priority_list = [
            "Protein kinase",
            #"Histidine kinase",
        ]
        i = 0
        while i < len(priority_list):
            for f in features:
                if f["name"] == priority_list[i]:
                    return seq[f["start"] - 1:f["end"]]
            i = i + 1
        return None
    

def get_aloop_seq(seq: str):
    if type(seq) != str:
        return None
    else:
        aloop_start_motifs = ["DFG", "DLG"]
        aloop_end_motifs = ["APE", "SPE"]
        motifs_pairs = [(start, end) for start in aloop_start_motifs for end in aloop_end_motifs]
        substrings = []
        for s in motifs_pairs:
            candidate_substr = extract_all_between(seq, s[0], s[1])
            if len(candidate_substr) > 0:
                candidate_substr = min(candidate_substr, key=lambda x: abs(len(x) - 30))
                substrings.append(f"{s[0]}{candidate_substr}{s[1]}")
        if len(substrings) > 0:
            aloop = min(substrings, key=lambda x: abs(len(x) - 30))
            if len(aloop) < 100:
                return aloop
        else:
            return None


def extract_all_between(text, pattern_1, pattern_2):
    pattern = f"{re.escape(pattern_1)}(.*?){re.escape(pattern_2)}"
    matches = re.findall(pattern, text, re.DOTALL)
    return matches


def main():
    kinases_tsv_gz = sys.argv[1]

    df = pd.read_csv(kinases_tsv_gz, sep="\t")

    df['Kinase domain'] = df.apply(
        lambda row: get_kinase_domain_seq(row['Sequence'], row['Domain [FT]']),
        axis=1
    )

    df['A-loop'] = df.apply(
        lambda row: get_aloop_seq(row['Kinase domain']),
        axis=1
    )

    df = df.replace({None: np.nan})
    
    df = df.drop(columns=["Entry Name", "Gene Names", "Domain [FT]"])

    print(df.to_csv(sep="\t", index=False, header=True))


if __name__ == "__main__":
    main()
