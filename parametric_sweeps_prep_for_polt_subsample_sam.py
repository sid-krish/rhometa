import os

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

def collect_results_sweep_1(rho, theta, sample_size, depth, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Load data into dataframe
    recom_est_results_dir = f"{os.getcwd()}/SimulationResults(subsample)/Recom_Est_Output/"

    col_names = ["rho_sim", "theta_sim", "sample_size_sim", "depth_sim", "seed_sim",
                 "mean", "std", "min", "5%", "25%", "50%", "75%", "95%", "max"]

    df_recom_est_resutls = pd.DataFrame(columns=col_names)

    for rho, theta, sample_size, depth, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_sample_size_{sample_size}_depth_{depth}_seed_{seed}_"

        if os.path.isfile(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv"):
            with open(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv", 'r') as results:
                results_unfiltered = results.readlines()[1].strip().split(',')
                results_final = results_unfiltered[2:]

            to_append = [rho, theta, sample_size, depth, seed] + results_final
            df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

        else:
            print(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv")
            results_final = 9 * [np.nan]
            to_append = [rho, theta, sample_size, depth, seed] + results_final
            df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

    return df_recom_est_resutls


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.01, 0.025, 0.05, 0.075, 0.1]
    theta_sweep_1 = [0.01]
    sample_size_sweep_1 = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    depth_sweep_1 = [.5, 1, 2, 4, 8, 10]
    seed_sweep_1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    recom_tract_len = 500

    collected_results_sweep_1_df = collect_results_sweep_1(rho_sweep_1, theta_sweep_1, sample_size_sweep_1,
                                                           depth_sweep_1, seed_sweep_1)

    # process and export df for plotting
    # Since fastsimbac does 2n*r only scaling by tract len is needed

    collected_results_sweep_1_df["scaled_rho_sim"] = collected_results_sweep_1_df["rho_sim"].apply(
        lambda x: x * recom_tract_len)

    reorder_cols = ['rho_sim', 'scaled_rho_sim', 'theta_sim', 'sample_size_sim', 'depth_sim', 'seed_sim',
                    'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']

    collected_results_sweep_1_df = collected_results_sweep_1_df.reindex(columns=reorder_cols)

    collected_results_sweep_1_df.to_csv("collected_results_sweep_1_subsample.csv", index=None)

    collected_results_sweep_1_df = collected_results_sweep_1_df.astype('float64')

    # Plot results
    sns.set_theme(style="whitegrid")

    ax = sns.catplot(data=collected_results_sweep_1_df, x="scaled_rho_sim", y="mean", hue="sample_size_sim",
                     col="depth_sim", col_wrap=3, sharex=True, sharey=True, kind="box", palette="YlGnBu")

    ax.set(ylim=(0, 50))

    ax.set(yticks=([5, 12.5, 25, 37.5, 50]))
    # ax.map(plt.axhline, y=5, ls='dotted', color='b', linewidth=1)
    # ax.map(plt.axhline, y=12.5, ls='--', color='black', linewidth=1)
    # ax.map(plt.axhline, y=25, ls='--', color='black', linewidth=1)
    # ax.map(plt.axhline, y=37.5, ls='--', color='black', linewidth=1)
    # ax.map(plt.axhline, y=50, ls='--', color='black', linewidth=1)

    ax.set(xlabel='Simulated Rho', ylabel='Estimated Rho (mean)')

    # ax.set_title("MPRR Simulated (scaled_rho_sim) vs Estimated Rho (mean)")

    # ax.figure.savefig("MPRR_results.png", dpi=500)

    ax.savefig("MPRR_subsample_results_mixed_sam.png", dpi=500)

    #### deviation plot

    collected_results_sweep_1_df["deviation"] = (collected_results_sweep_1_df["mean"] -
                                                 collected_results_sweep_1_df["scaled_rho_sim"]) / collected_results_sweep_1_df[
                                                    "scaled_rho_sim"]

    ax2 = sns.catplot(data=collected_results_sweep_1_df, x="scaled_rho_sim", y="deviation", hue="sample_size_sim",
                     col="depth_sim", col_wrap=3, sharex=True, sharey=True, kind="box", palette="RdYlGn")

    ax2.set(ylim=(-1,1))

    ax2.map(plt.axhline, y=0, ls='--', color='g', linewidth=2)

    ax2.savefig("MPRR_subsample_deviation_mixed_sam.png", dpi=500)
