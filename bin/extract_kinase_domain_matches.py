#!/usr/bin/env python3

import sys
from Bio import SeqIO


"""
domtblout_file = "input/kinase_domains.hmmsearch"
fasta_file = "input/kinases.fa"
output_file = "kinase_metadata/kinase_domains_hmm_matches.fa"
"""


def main():

    domtblout_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    # Step 1: Parse domtblout to get intervals
    intervals = {}
    with open(domtblout_file, "r") as f:
        header_count = 0
        in_a_header = False
        for line in f:
            if header_count < 2 and not line.startswith("#"):  # Skip header/comment lines
                in_a_header = False
                fields = line.split()
                target = fields[0]
                start = int(fields[17])  # ali coord from
                end = int(fields[18])    # ali coord to
                # Store intervals (note: start is 1-based, end is inclusive)
                if target not in intervals:
                    intervals[target] = []
                intervals[target].append((start, end))
            elif header_count >= 2:
                break
            else:
                if in_a_header == False:
                    header_count = header_count + 1
                    in_a_header = True

    # Step 2: Load FASTA file and extract subsequences
    with open(output_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            target = record.id
            if target in intervals:
                sequence = str(record.seq)
                for start, end in intervals[target]:
                    # Adjust for 0-based indexing (subtract 1 from start)
                    subsequence = sequence[start-1:end]
                    # Write to output in FASTA format
                    out.write(f">{target}_{start}-{end}\n{subsequence}\n")


if __name__ == "__main__":
    main()
