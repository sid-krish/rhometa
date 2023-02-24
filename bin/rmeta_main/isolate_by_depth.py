#!/usr/bin/env python


def main(pairwise_table, depth):
    df = pairwise_table

    # Isolate slice of data based on depth
    df = df.loc[df.sum(axis=1) == depth]

    return df
