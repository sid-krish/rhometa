#!/bin/bash
#PBS -M 11849395@uts.edu.au -m abe
#PBS -l select=1:ncpus=2:mem=8GB -l walltime=120:00:00

cd $PBS_O_WORKDIR

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate

nextflow run align_reads.nf -qs 50 -resume