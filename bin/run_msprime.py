#!/usr/bin/env python

import msprime
import numpy as np

sample_size = 50
Ne = 1
length = 100000
mutation_rate = 0.01
recombination_rate = 0.000025
random_seed = 123


tree_sequence = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, mutation_rate=mutation_rate,
                                 recombination_rate=recombination_rate, random_seed=random_seed)

print(f"rho (4*N0*r*l) = {4*Ne*recombination_rate*length}")

with open("cleanTrees.txt", 'w') as f:
    for tree in tree_sequence.trees():
        newick_encoded_tree = tree.newick() # get newick using tskit library. Tree is a tskit object
        # Have to do some further processing for seq-gen
        # https://darencard.net/blog/2019-10-22-msprime-tutorial/
        tree_length = -int(np.round(tree.get_interval()[0]))+int(np.round(tree.get_interval()[1]))

        f.write(f"[{tree_length}]{newick_encoded_tree}")
        f.write('\n')
