#!/usr/bin/env python
import msprime
import tskit

ts = msprime.sim_ancestry(samples=50, ploidy=1, sequence_length=100000,
                          gene_conversion_rate=0.005, gene_conversion_tract_length=1000, random_seed=1234)  # rho 10, final result should be scaled by 2
mts = msprime.sim_mutations(ts, rate=0.01, random_seed=1234) # Final result should be scaled by 2

ref_seq = tskit.random_nucleotides(length=100000, seed=1234)
mts.write_fasta(file_or_path="msp_rho_test10.fa",
                reference_sequence=ref_seq, wrap_width=0)
