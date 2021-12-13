import os
import subprocess

import numpy as np
from tqdm import tqdm

"""
TODO: 2 Parametric sweeps + Replicates

####
Sweep 1: Recombination rate simulations
####

Rho: 0.1, 0.5, 2.5, 25, 50
Theta: 0.01 (fixed like PIIM)
Genome_size: 100,000; 250,000; 500,000; 750,000; 1,000,000
mean_cov: .5, 1, 2.5, 5, 10
Samples: 10 (currently fixed)

####
Sweep 2: Mutation rate simulations
####

Theta: 0.001, 0.005, 0.01, 0.05, 0.1
Rho: 0.01 (fixed)
Genome_size: 100,000; 250,000; 500,000; 750,000; 1,000,000
mean_cov: .5, 1, 2.5, 5, 10
Samples: 10 (currently fixed)

####
3 Replicates with different seed values
####

Seed: 123, 456, 789

"""


def sweep_1_simulation(rho, theta, sample_size, depth, genome_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, genome_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 6)

    # Generate simulated datasets
    for rho, theta, sample_size, depth,genome_size, seed in tqdm(sweep_1_combinations):
        subprocess.run(["nextflow", "run", "sim_gen.nf",
                        "--rho_rates", f"{rho}",
                        "--mutation_rate", f"{theta}",
                        "--sample_sizes", f"{sample_size}",
                        "--fold_cov", f"{depth}",
                        "--genome_sizes", f"{genome_size}",
                        "--seed", f"{seed}",
                        "--prepend_filename",
                        f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_genome_size_{genome_size}_seed_{seed}_"])

    return None


def sweep_1_recom_est(rho, theta, sample_size, depth, genome_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, genome_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 6)

    # Process simulated datasets
    simulation_results_dir = f"{os.getcwd()}/Sim_Gen_Output/"
    for rho, theta, sample_size, depth,genome_size, seed in tqdm(sweep_1_combinations):
        prepended_filename = f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_genome_size_{genome_size}_seed_{seed}_"
        subprocess.run(["nextflow", "run", "recom_est.nf",
                        "--bam_file",
                        f"{simulation_results_dir}{prepended_filename}Aligned.bam",
                        "--reference_genome",
                        f"{simulation_results_dir}{prepended_filename}firstGenome.fa",
                        "--prepend_filename",
                        f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_genome_size_{genome_size}_seed_{seed}_",
                        "--subsample_bam"])

    return None


def sweep_2_simulation(rho, theta, sample_size, depth, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Generate simulated datasets
    for rho, theta, sample_size, depth, seed in sweep_1_combinations:
        subprocess.run(["nextflow", "run", "sim_gen.nf",
                        "--rho_rates", f"{rho}",
                        "--mutation_rate", f"{theta}",
                        "--sample_sizes", f"{sample_size}",
                        "--fold_cov", f"{depth}",
                        "--seed", f"{seed}",
                        "--prepend_filename",
                        f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_seed_{seed}_"])

    return None


def sweep_2_mut_est(rho, theta, sample_size, depth, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Process simulated datasets
    simulation_results_dir = f"{os.getcwd()}/Sim_Gen_Output/"
    for rho, theta, sample_size, depth, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_seed_{seed}_"
        subprocess.run(["nextflow", "run", "theta_est.nf",
                        "--bam_file",
                        f"{simulation_results_dir}{prepended_filename}Aligned.bam",
                        "--reference_genome",
                        f"{simulation_results_dir}{prepended_filename}firstGenome.fa",
                        "--prepend_filename",
                        f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_seed_{seed}_"])

    return None


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.01, 0.025, 0.05, 0.075, 0.1]
    theta_sweep_1 = [0.01]
    sample_size_sweep_1 = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    depth_sweep_1 = [.5, 1, 2, 4, 8, 10]
    genome_size_sweep_1 = [25000]
    seed_sweep_1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    sweep_1_simulation(rho_sweep_1, theta_sweep_1, sample_size_sweep_1, depth_sweep_1,genome_size_sweep_1, seed_sweep_1)
    sweep_1_recom_est(rho_sweep_1, theta_sweep_1, sample_size_sweep_1, depth_sweep_1,genome_size_sweep_1, seed_sweep_1)

    # Sweep 2: Mutation rate estimation
    # theta_sweep_2 = [0.001, 0.005, 0.01, 0.05, 0.1]
    # rho_sweep_2 = [0.01]
    # sample_size_sweep_2 = [100000, 250000, 500000, 750000, 1000000]
    # depth_sweep_2 = [.5, 1, 2.5, 5, 10]
    # seed_sweep_2 = [123, 456, 789]

    # sweep_2_simulation(rho_sweep_2, theta_sweep_2, sample_size_sweep_2, depth_sweep_2, seed_sweep_2)
    # sweep_2_mut_est(rho_sweep_2, theta_sweep_2, sample_size_sweep_2, depth_sweep_2, seed_sweep_2)
