#!/usr/bin/env python

import sys

import pandas as pd

pileup_table = sys.argv[1]

df = pd.read_csv(pileup_table, index_col="RefPos_0-based")

# cols_to_drop = ["Num_bases_aligned"]

# df = df.drop(columns=cols_to_drop)

df = df[df["Allele_counts"] > 1]

df.to_csv("variants_from_pileup.csv")