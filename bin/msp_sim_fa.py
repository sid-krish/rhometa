#!/usr/bin/env python
import sys
import msprime
import tskit

samples_input = int(sys.argv[1])
sequence_length_input = int(sys.argv[2])
gene_conversion_rate_input = float(sys.argv[3])
gene_conversion_tract_length_input = int(sys.argv[4])
random_seed_input = int(sys.argv[5])
mutation_rate = float(sys.argv[6])

ts = msprime.sim_ancestry(  # Final result is recombination rate * 2 * N_e * tract length
    samples=samples_input,
    ploidy=1,  # Fixed. Only interested in haploid samples
    sequence_length=sequence_length_input,
    gene_conversion_rate=gene_conversion_rate_input,
    gene_conversion_tract_length=gene_conversion_tract_length_input,
    random_seed=random_seed_input,
)  # population_size "If not specified, defaults to 1"

mts = msprime.sim_mutations(  # Final result is mutataion rate * N_e * 2
    ts, rate=mutation_rate, random_seed=random_seed_input
)


ref_seq = tskit.random_nucleotides(length=sequence_length_input, seed=random_seed_input)

mts.write_fasta(file_or_path="msp_out.fa", reference_sequence=ref_seq, wrap_width=0)
