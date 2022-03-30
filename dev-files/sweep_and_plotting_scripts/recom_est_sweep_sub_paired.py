import os
import subprocess

import numpy as np
from tqdm import tqdm


def recom_est(rho, theta, sample_size, depth, genome_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, genome_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    combinations = mesh_grid.T.reshape(-1, 6)

    # Process simulated datasets
    simulation_results_dir = f"{os.getcwd()}/Sim_Gen_Output/"
    for rho, theta, sample_size, depth, genome_size, seed in tqdm(combinations):
        prepended_filename = f"rho_{rho}_theta_{theta}_sample_size_{int(sample_size)}_depth_{int(depth)}_genome_size_{int(genome_size)}_seed_{int(seed)}_"
        subprocess.run(["nextflow", "run", "recom_est.nf",
                        "--bam_file",
                        f"{simulation_results_dir}{prepended_filename}Aligned.bam",
                        "--reference_genome",
                        f"{simulation_results_dir}{prepended_filename}firstGenome.fa",
                        "--subsample_bam"])

    return None


if __name__ == '__main__':
    rho_sweep_1 = [0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.005, 0.015, 0.025, 0.035, 0.045]
    theta_sweep_1 = [0.01]
    sample_size_sweep_1 = [10, 20, 40, 80, 160, 200]
    depth_sweep_1 = [1, 2, 4, 8]
    genome_size_sweep_1 = [100000]
    seed_sweep_1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    recom_est(rho_sweep_1, theta_sweep_1, sample_size_sweep_1, depth_sweep_1, genome_size_sweep_1, seed_sweep_1)
