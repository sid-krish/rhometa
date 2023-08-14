#!/bin/bash
#PBS -M 11849395@uts.edu.au -m abe
#PBS -l select=1:ncpus=2:mem=12g -l walltime=100:00:00

. ~/miniconda3/etc/profile.d/conda.sh
conda activate nextflow

cd $PBS_O_WORKDIR

nextflow run viral_analysis.nf --contigs $CTG --out_dir $OUT