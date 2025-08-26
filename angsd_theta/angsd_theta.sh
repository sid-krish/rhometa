# Need to think about maths behind these steps so that it can be documented. Which parts are ML and which parts are Bayesian?

# Commands follow guide from https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests

#First estimate the site allele frequency likelihood
angsd -i /home/sid/Github/rhometa/testing/recom_0.005_tract_1000_mutation_0.005_sample_size_20_depth_4_genome_size_50000_seed_123_final.bam -doSaf 1 -anc recom_0.005_tract_1000_mutation_0.005_sample_size_20_depth_4_genome_size_50000_seed_123_final.fa -GL 1 -P 4 -out out #output is log scale

#For folded " If you do not have the ancestral state you can simply use the assembly you have mapped agains, but remember to add -fold 1 in the 'realSFS' and 'realSFS sf2theta' step."

# Obtain the maximum likelihood estimate of the SFS using the realSFS program
realSFS out.saf.idx -P 4 -fold 1 > out.sfs

# Calculate the thetas for each site
realSFS saf2theta out.saf.idx -outname out -sfs out.sfs -fold 1 
# Estimate Tajimas D and other statistics
thetaStat do_stat out.thetas.idx # no longer log scale. Final vals. pestPG file are the sum of the per site estimates for a region -> chromosome wide vals

# Sliding wnindow analysis. Not needed. Testing
# thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz