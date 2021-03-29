#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys


def fourNonZeroValsFilter(df):
    """
    Isolate rows where 4 different pairs of bases are detected
    :return: filter_1_df
    """

    df = df.replace(to_replace=0, value=np.nan)  # replace 0 in df with nan for df.count step

    df["NonZeroVals"] = df.count(axis="columns")  # Count cols (values in a row) that are not nan

    # Get rows where value in "NonZeroVals" col is 4 (the rows that fail the four gamete test)
    filtered_df = df.loc[df["NonZeroVals"] == 4]

    filtered_df = filtered_df.drop(columns="NonZeroVals") # Must be dropped for next step
    return filtered_df


def find_seq_error(row):
    rowVals = row.values
    if 1 in rowVals:
        return 'y'

    return 'n'


def seqErrorFilter(df):
    # if val 1 in row those rows need to be removed, as those values likely represent sequence errors
    df["seq_error"] = df.apply(find_seq_error, axis=1)  # axis = 1, apply function to row
    filtered_df = df.loc[df["seq_error"] == 'n'] # isolate rows that don't have 1
    filtered_df = filtered_df.drop(columns="seq_error")  # No longer needed

    return filtered_df


def getBiallelicSets(row):
    entry = tuple(zip(row.index,row.values))

    #Step 1: Represent all pairs with values as one merged string
    pairsWithVals = []

    for pair, value in entry:
        if str(value) != "nan":
            pairsWithVals.append((pair))

    mergedPairs = str(''.join(pairsWithVals))

    #Step 2: Get count of unique bases in site 1 and site 2
    site1s = mergedPairs[0::2] # Site 1 bases
    unqSite1Vals = len(set(site1s)) #Get count of unique values in site 1 bases

    site2s = mergedPairs[1::2]
    unqSite2Vals = len(set(site2s))

    #Step 3: If only 2 unique bases in site 1 and site 2 return confirmation
    if unqSite1Vals == 2 and unqSite2Vals == 2:
        return 'y'

    return 'n'


def filter_bi_allelic(df):
    #Ensure site 1 and site 2 bases are biallelic
    df["Bi-allelic_sets"] = df.apply(getBiallelicSets, axis = 1) # axis = 1, apply function to row

    filtered_df = df.loc[df["Bi-allelic_sets"] == 'y']

    filtered_df = filtered_df.drop(columns="Bi-allelic_sets")  # No longer needed
    return filtered_df


def export_final_df(final_df):

    final_df.to_csv("FGT_results.csv")

    return None


if __name__ == '__main__':

    fgtTable = sys.argv[1]

    df = pd.read_csv(fgtTable, index_col="RefPos")

    df_filter_1 = fourNonZeroValsFilter(df)
    df_filter_2 = seqErrorFilter(df_filter_1)
    df_filter_3 = filter_bi_allelic(df_filter_2)

    export_final_df(df_filter_3)
