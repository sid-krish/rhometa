#!/bin/bash
#PBS -M 11849395@uts.edu.au -m abe
#PBS -l select=1:ncpus=12:mem=20GB:host=i3node01 -l walltime=120:00:00 -q i3q

cd $PBS_O_WORKDIR

nextflow run sub_theta_align_est.nf -qs 50 -resume