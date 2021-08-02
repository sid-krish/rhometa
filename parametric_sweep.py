import subprocess
import numpy as np

"""
TODO: 2 Parametric sweeps + Replicates

####
Sweep 1: Recombination rate simulations
####

Rho: 10, 30, 50, 60, 90
Theta: 0.01 (fixed like PIIM)
Genome_size: 100,000; 250,000; 500,000; 750,000; 1,000,000
Depth: 5, 25, 50, 75, 100

####
Sweep 2: Mutation rate simulations
####

Theta: 0.1, 0.3, 0.5, 0.6, 0.9
Genome_size: 100,000; 250,000; 500,000; 750,000; 1,000,000
Depth: 5, 25, 50, 75, 100

####
3 Replicates with different seed values
####

Seed: 123, 456, 789

"""


def sweep_1(rho, theta, genome_size, depth, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, genome_size, depth, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Generate simulated datasets
    for rho, theta, genome_size, depth, seed in sweep_1_combinations:
        subprocess.run(["nextflow", "run", "sim_gen.nf",
                        "--rho_rates", f"{rho}",
                        "--mutation_rate", f"{theta}",
                        "--genome_sizes", f"{genome_size}",
                        "--fold_cov", f"{depth}",
                        "--seed", f"{seed}",
                        "--prepend_filename",
                        f"rho_{rho}_theta_{theta}_genome_size_{genome_size}_depth_{depth}_seed_{seed}_"])

    # for rho estimator get performance stats
    return None


if __name__ == '__main__':

    # Sweep 1: Recombination rate simulations
    rho = [10, 30, 50, 60, 90]
    theta = [0.01]
    genome_size = [100000, 250000, 500000, 750000, 1000000]
    depth = [5, 25, 50, 75, 100]
    seed = [123, 456, 789]

    sweep_1(rho, theta, genome_size, depth, seed)
