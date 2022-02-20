import os

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

def collect_results_sweep_1(rho, theta, sample_size, depth, genome_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, sample_size, depth, genome_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 6)

    # Load data into dataframe
    theta_est_results_dir = f"/Users/Sid/Documents/Github/rhometa/Misc/msp_theta_sweep/Theta_Est_Output/"

    col_names = ["rho_sim", "theta_sim", "sample_size_sim", "depth_sim", "genome_size_sim", "seed_sim",
                 "mean_depth", "theta_per_site_at_mean_depth", "median_depth", "theta_per_site_at_median_depth"]

    df_theta_est_results = pd.DataFrame(columns=col_names)

    for rho, theta, sample_size, depth, genome_size, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_sample_size_{int(sample_size)}_depth_{int(depth)}_genome_size_{int(genome_size)}_seed_{int(seed)}_final_"

        if os.path.isfile(f"{theta_est_results_dir}{prepended_filename}Theta_estimate_stats.csv"):
            df = pd.read_csv(f"{theta_est_results_dir}{prepended_filename}Theta_estimate_stats.csv", header=None)
            df = df.T

            results_final = df.iloc[1].to_list()

            to_append = [rho, theta, sample_size, depth, genome_size, seed] + results_final
            df_theta_est_results.loc[len(df_theta_est_results)] = to_append

        else:
            print(f"{theta_est_results_dir}{prepended_filename}Theta_estimate_stats.csv")
            results_final = 4 * [np.nan]
            to_append = [rho, theta, sample_size, depth, genome_size, seed] + results_final
            df_theta_est_results.loc[len(df_theta_est_results)] = to_append

    return df_theta_est_results


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.0]
    theta_sweep_1 = [0.0, 0.05, 0.15, 0.25, 0.35] # unscaled u values. theta = 2 . p . N_e . u
    sample_size_sweep_1 = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
    depth_sweep_1 = [1, 4, 8, 16]
    genome_size_sweep_1 = [100000]
    seed_sweep_1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

    recom_tract_len = 1000

    collected_results_sweep_1_df = collect_results_sweep_1(rho_sweep_1, theta_sweep_1, sample_size_sweep_1, depth_sweep_1,genome_size_sweep_1, seed_sweep_1)

    # process and export df for plotting
    # Since fastsimbac does 2n*r only scaling by tract len is needed

    collected_results_sweep_1_df["scaled_theta_sim"] = collected_results_sweep_1_df["theta_sim"].apply(
        lambda x: x * 2)

    reorder_cols = ["rho_sim", "theta_sim", "scaled_theta_sim", "sample_size_sim", "depth_sim", "genome_size_sim", "seed_sim",
                     "mean_depth", "theta_per_site_at_mean_depth", "median_depth", "theta_per_site_at_median_depth"]

    collected_results_sweep_1_df = collected_results_sweep_1_df.reindex(columns=reorder_cols)

    collected_results_sweep_1_df.to_csv("collected_theta_results_subsample.csv", index=None)

    collected_results_sweep_1_df = collected_results_sweep_1_df.astype('float64')

    collected_results_sweep_1_df.rename(columns = {'sample_size_sim':'genomes', 'depth_sim':'fold_coverage'}, inplace = True) # Rename for plot

    # Plot results
    sns.set_theme(style="white")

    # x-axis variable is treated as categorical

    ax = sns.catplot(data=collected_results_sweep_1_df, x="scaled_theta_sim", y='theta_per_site_at_median_depth', hue="genomes",
                     col="fold_coverage", col_wrap=2, sharex=True, sharey=True, palette="mako_r", kind="box", linewidth=0.2)

    # ax.set(ylim=(0, 50), xlabel="Simulated \u03C1", ylabel="Estimated \u03C1 (median)")

    # ax.set(yticks=([1.0, 5.0, 15.0, 25.0, 35.0, 45.0]))
    # ax.set(yticks=(range(0,45)))

    # ax.map(plt.axhline, y=0.1, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=0.25, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=0.5, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=0.75, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=1, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=5, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=15, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=25, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=35, ls='dotted', color='g', linewidth=1)
    # ax.map(plt.axhline, y=45, ls='dotted', color='g', linewidth=1)

    # ax.fig.suptitle("rhometa Simulated (scaled_rho_sim) vs Estimated Rho (median)", y=1.05)

    # ax.figure.savefig("rhometa_results.png", dpi=500)

    ax.savefig("rhometa_theta_est_results.png", dpi=500)

    #### deviation plot

    collected_results_sweep_1_df["deviation"] = (collected_results_sweep_1_df['theta_per_site_at_median_depth'] -
                                                 collected_results_sweep_1_df["scaled_theta_sim"]) / collected_results_sweep_1_df[
                                                    "scaled_theta_sim"]
    
    ax2 = sns.catplot(data=collected_results_sweep_1_df, x="scaled_theta_sim", y="deviation", hue="genomes",
                     col="fold_coverage", col_wrap=2, sharex=True, sharey=True, kind="box", palette="RdYlGn", linewidth=0.2)
    
    ax2.set(ylim=(-1,1), xlabel="Simulated \u03C1", ylabel="Deviation")
    
    ax2.map(plt.axhline, y=0, ls='--', color='g', linewidth=2)
    
    ax2.savefig("rhometa_theta_est_deviation_results.png", dpi=500)
