#!/usr/bin/env python
import sys

import pandas as pd


def getBiallelicSets(row):
    # row.drop("d_ij", inplace=True)
    entry = tuple(zip(row.index, row.values))

    # Step 1: Represent all pairs with values as one merged string
    pairsWithVals = []

    for pair, value in entry:
        if int(value) != 0:
            pairsWithVals.append((pair))

    mergedPairs = str(''.join(pairsWithVals))
    # Step 2: Get count of unique bases in site 1 and site 2
    site1s = mergedPairs[0::2]  # Site 1 bases
    unqSite1Vals = len(set(site1s))  # Get count of unique values in site 1 bases

    site2s = mergedPairs[1::2]
    unqSite2Vals = len(set(site2s))

    # Step 3: If only 2 unique bases in site 1 and site 2 return confirmation
    if unqSite1Vals == 2 and unqSite2Vals == 2:
        return 'y'

    return 'n'


def filter_bi_allelic(df):
    # Ensure site 1 and site 2 bases are biallelic
    df["Bi-allelic_sets"] = df.apply(getBiallelicSets, axis=1)  # axis = 1, apply function to row

    filtered_df = df.loc[df["Bi-allelic_sets"] == 'y']

    filtered_df = filtered_df.drop(columns="Bi-allelic_sets")  # No longer needed
    return filtered_df


def export_final_df(final_df, depth):
    final_df.to_csv(f"pairwise_biallelic_table_depth_{depth}.csv")

    return None


if __name__ == '__main__':

    pairwise_table = sys.argv[1]

    depth = int(sys.argv[2])

    # pairwise_table = "pairwise_table.csv"

    df = pd.read_csv(pairwise_table, index_col="RefPos_0-based")

    # Isolate slice of data based on depth
    df = df[df.sum(axis=1) == depth]

    df_filtered = filter_bi_allelic(df)

    if not df_filtered.empty:
        export_final_df(df_filtered, depth)
