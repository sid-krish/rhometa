import os

import numpy as np
import pandas as pd
import seaborn as sns


def collect_results_sweep_1(rho, theta, genome_size, depth, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, genome_size, depth, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Load data into dataframe
    recom_est_results_dir = f"{os.getcwd()}/SimulationResults/Recom_Est_Output/"

    col_names = ["rho_sim", "theta_sim", "genome_size_sim", "depth_sim", "seed_sim",
                 "mean","std","min","5%","25%","50%","75%","95%","max"]

    df_recom_est_resutls = pd.DataFrame(columns=col_names)

    for rho, theta, genome_size, depth, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_genome_size_{genome_size}_depth_{depth}_seed_{seed}_"
        with open(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv", 'r') as results:
            results_unfiltered = results.readlines()[1].strip().split(',')
            results_final = results_unfiltered[2:]

        to_append = [rho, theta, genome_size, depth, seed] + results_final
        df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

    return df_recom_est_resutls


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation plot
    rho_sweep_1 = [0.1, 0.5, 2.5, 25, 50]
    theta_sweep_1 = [0.01]
    genome_size_sweep_1 = [100000, 250000, 500000, 750000, 1000000]
    depth_sweep_1 = [.5, 1, 2.5, 5, 10]
    seed_sweep_1 = [123, 456, 789]

    collected_results_sweep_1_df = collect_results_sweep_1(rho_sweep_1, theta_sweep_1, genome_size_sweep_1, depth_sweep_1, seed_sweep_1)

    # plot_results_sweep_1()
