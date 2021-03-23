#!/usr/bin/env python
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO


def get_size_of_largest_genome(input_fasta):
    records = SeqIO.parse(input_fasta, "fasta")  # gets exhausted once accessed

    genome_sizes = [len(record.seq) for record in records]

    genome_sizes.sort()

    largest = genome_sizes.pop()

    return int(largest)


def convert_fasta_to_2D_matrix(input_fasta, size_of_largest_genome):
    df = pd.DataFrame(index=range(size_of_largest_genome))

    records = SeqIO.parse(input_fasta, "fasta")  # gets exhausted once accessed

    for record in records:
        genome_array = np.array(record.seq, dtype=np.str_)
        df[record.description] = genome_array

    return df


def get_alleles(df_fasta):
    """
    This method can be parallelized using dask
    """

    df_fasta["alleles"] = df_fasta.apply(lambda x: pd.unique(x), axis=1)  # axis = 1, apply function to row

    return df_fasta


def get_allele_counts(df_fasta):
    """
    This method can be parallelized using dask
    """

    df_fasta["allele_counts"] = df_fasta.apply(lambda x: len(x["alleles"]), axis=1)  # axis = 1, apply function to roww

    return df_fasta


def main():
    # if __name__ == '__main__':

    input_fasta = sys.argv[1]
    # input_fasta = "../Output/s_123_m_0.01_p_20/s_123_m_0.01_p_20_reformated.fa"

    size_of_largest_genome = get_size_of_largest_genome(input_fasta)

    # switching to dask here can alleviate memory issues
    df_fasta = convert_fasta_to_2D_matrix(input_fasta, size_of_largest_genome)

    df_fasta = get_alleles(df_fasta)

    df_fasta = get_allele_counts(df_fasta)

    df_fasta.index.name = "RefPos_0-based"  # rename index

    # df_fasta.to_csv("ref_variants_in_fasta.csv")

    df_fasta = df_fasta[df_fasta["allele_counts"] > 1]

    drop_cols = df_fasta.columns.values[:-2]  # all cols except last 2

    df_fasta = df_fasta.drop(columns=drop_cols)

    df_fasta.to_csv("variants_in_fasta.csv")


if __name__ == '__main__':
    main()

"""
its possible to do this row by row, with significantly less memory use by avoiding pandas
Not sure how speed will be impacted though

pseudo code:

max_size = get_size_of_largest_genome()

variants = []

for refPos in range(max_size):
    alleles = []
    for genome in fasta:
        alleles = genome.base[i].append()
    
    unique_alleles = alleles.unique()

    alleles_count = len(unique_alleles)

    variants.append([i,unique_alleles,alleles_count])

"""
