#!/bin/bash
#PBS -M 11849395@uts.edu.au -m abe
#PBS -l select=1:ncpus=12:mem=30GB -l walltime=120:00:00 -q i3q

cd $PBS_O_WORKDIR

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate

nextflow run combined_align_est.nf -qs 50 -resume