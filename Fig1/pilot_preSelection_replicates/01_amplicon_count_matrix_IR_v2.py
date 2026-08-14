# Notes: R2 should be slightly more reliable as it's the read where the variable sequence is closer to the start of the read.

import re
import numpy as np
import pandas as pd
import pyfastx


# ------------- Input variables ------------- #

flank_5prime = "TCTAGT"
flank_3prime = "TATTCG"
varseq_length = 9
possible_seqs_path = "inputs/variants.csv"
samples_path = "inputs/samples_trimmed.tsv"


# ------------- Auxiliary functions ------------- #

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


# ------------- Setup ------------- #

search_seq = f"{flank_5prime}.{{{varseq_length}}}{flank_3prime}"
search_seq_reverse = f"{reverse_complement(flank_3prime)}.{{{varseq_length}}}{reverse_complement(flank_5prime)}"

# Pre-compile regex patterns once
pattern_fwd = re.compile(search_seq)
pattern_rev = re.compile(search_seq_reverse)

# Slice indices for extracting the variable sequence from a forward match
# e.g. "TCTAGT" is 6bp, so variable region starts at index 6
var_start = len(flank_5prime)
var_end = var_start + varseq_length

samples = pd.read_csv(samples_path, sep="\t", header=None)
variants = pd.read_csv(possible_seqs_path, header=None)

# Use a set for O(1) variant lookups instead of repeated .tolist() calls
variant_set = set(variants[0])
sample_names = samples[0].tolist()

# Use nested dicts for counts instead of pandas .loc (much faster for incremental updates)
# Structure: {sample_name: {variant: count}}
count_dict = {s: {v: 0 for v in variant_set} for s in sample_names}

# Accumulate QC rows in a list; concat once at the end
qc_rows = []

# ------------- Main loop ------------- #

for _, row in samples.iterrows():
    sample_name, r1_path, r2_path = row[0], row[1], row[2]

    no_hit_varseq = 0
    flanking_seq_errors = 0
    var_seq_mismatch = 0

    sample_counts = count_dict[sample_name]  # local reference avoids repeated dict lookups

    for seq1, seq2 in zip(pyfastx.Fastx(r1_path), pyfastx.Fastx(r2_path)):
        if seq1[0] != seq2[0]:
            print(f"Fastq file indexes got off-track from one another at sample: {sample_name}")
            print(f"At R1 index: {seq1[0]} and R2 index: {seq2[0]}")
            break

        r1_match = pattern_fwd.search(seq1[1])
        r2_match = pattern_rev.search(seq2[1])

        if r1_match and r2_match:
            r1_val = r1_match.group(0)
            # Compare forward match against reverse-complemented reverse match
            if r1_val == reverse_complement(r2_match.group(0)):
                final_val = r1_val[var_start:var_end]
                if final_val in variant_set:
                    sample_counts[final_val] += 1
                else:
                    no_hit_varseq += 1
            else:
                var_seq_mismatch += 1
        else:
            flanking_seq_errors += 1

    mapped = sum(sample_counts.values())
    qc_rows.append({
        "Sample": sample_name,
        "Mapped_reads": mapped,
        "No_flanking_sequence_hit": flanking_seq_errors,
        "R1_R2_var_seqs_mismatched": var_seq_mismatch,
        "var_seq_no_ref_match": no_hit_varseq,
    })

# ------------- Build outputs ------------- #

# Build count matrix from dict (single DataFrame construction, no repeated .loc writes)
count_matrix = pd.DataFrame(count_dict, index=sorted(variant_set)).reindex(columns=sample_names)

# Build QC table from accumulated rows (single concat)
qc = pd.DataFrame(qc_rows, columns=["Sample", "Mapped_reads", "No_flanking_sequence_hit",
                                      "R1_R2_var_seqs_mismatched", "var_seq_no_ref_match"])

count_matrix.to_csv("count_matrix.csv")
qc.to_csv("metrics.csv", index=False)
