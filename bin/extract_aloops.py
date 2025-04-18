#!/usr/bin/env python3

import sys
import re
from Bio import SeqIO


def extract_all_between(text, pattern_1, pattern_2):
    pattern = f"{re.escape(pattern_1)}(.*?){re.escape(pattern_2)}"
    matches = re.findall(pattern, text, re.DOTALL)
    return matches


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


def main():
    fasta_file = sys.argv[1]

    # Load FASTA file and look for a-loops heuristically
    for record in SeqIO.parse(fasta_file, "fasta"):
        target = record.id
        sequence = str(record.seq)
        aloop = get_aloop_seq(sequence)
        print(f">{target}")
        if aloop is None:
            print("none")
        else:
            print(aloop)


if __name__ == "__main__":
    main()
