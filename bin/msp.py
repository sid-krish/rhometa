#!/usr/bin/env python
import msprime
import tskit

ts = msprime.sim_ancestry(samples=10, ploidy=1, sequence_length=10000, gene_conversion_rate=0.2, gene_conversion_tract_length=100, random_seed=1234) #rho 20
mts = msprime.sim_mutations(ts, rate=0.03, random_seed=1234)

# ref_seq=tskit.random_nucleotides(length=10000, seed=1234)
mts.write_fasta(file_or_path="msp_rho_test.fa")