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
    recom_est_results_dir = f"/Users/Sid/Documents/Github/rhometa/Misc/simulation_results/Recom_Est_Output_fb(old)/"

    col_names = ["rho_sim", "theta_sim", "sample_size_sim", "depth_sim", "genome_size_sim", "seed_sim",
                 "mean", "std", "min", "5%", "25%", "50%", "75%", "95%", "max"]

    df_recom_est_resutls = pd.DataFrame(columns=col_names)

    for rho, theta, sample_size, depth, genome_size, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_sample_size_{int(sample_size)}_depth_{int(depth)}_genome_size_{int(genome_size)}_seed_{int(seed)}_final_"

        if os.path.isfile(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv"):
            with open(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv", 'r') as results:
                results_unfiltered = results.readlines()[1].strip().split(',')
                results_final = results_unfiltered[2:]

            to_append = [rho, theta, sample_size, depth, genome_size, seed] + results_final
            df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

        else:
            print(f"{recom_est_results_dir}{prepended_filename}final_results_summary.csv")
            results_final = 9 * [np.nan]
            to_append = [rho, theta, sample_size, depth, genome_size, seed] + results_final
            df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

    return df_recom_est_resutls


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.001, 0.005, 0.015, 0.025, 0.035, 0.045]
    # rho_sweep_1 = [0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.005, 0.015, 0.025, 0.035, 0.045]
    theta_sweep_1 = [0.01]
    sample_size_sweep_1 = [10, 20, 40, 80, 160, 200]
    depth_sweep_1 = [1, 2, 4, 8]
    genome_size_sweep_1 = [100000]
    seed_sweep_1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    recom_tract_len = 1000

    collected_results_sweep_1_df = collect_results_sweep_1(rho_sweep_1, theta_sweep_1, sample_size_sweep_1, depth_sweep_1,genome_size_sweep_1, seed_sweep_1)

    # process and export df for plotting
    # Since fastsimbac does 2n*r only scaling by tract len is needed

    collected_results_sweep_1_df["scaled_rho_sim"] = collected_results_sweep_1_df["rho_sim"].apply(
        lambda x: x * recom_tract_len)

    reorder_cols = ['rho_sim', 'scaled_rho_sim', 'theta_sim', 'sample_size_sim', 'depth_sim', 'seed_sim',
                    'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']

    collected_results_sweep_1_df = collected_results_sweep_1_df.reindex(columns=reorder_cols)

    collected_results_sweep_1_df.to_csv("collected_results_sweep_1_subsample.csv", index=None)

    collected_results_sweep_1_df = collected_results_sweep_1_df.astype('float64')

    collected_results_sweep_1_df.rename(columns = {'sample_size_sim':'Genomes', 'depth_sim':'fold_coverage_scaler'}, inplace = True) # Rename for plot

    # Plot results
    sns.set_theme(style="white")

    # x-axis variable is treated as categorical

    ax = sns.catplot(data=collected_results_sweep_1_df, x="scaled_rho_sim", y="mean", hue="Genomes",
                     col="fold_coverage_scaler", col_wrap=2, sharex=True, sharey=True, palette="Accent", kind="box")

    # ax.set(ylim=(0, 50), xlabel="Simulated \u03C1", ylabel="Estimated \u03C1 (mean)")

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

    ax.set(ylabel="Estimated \u03C1 (mean)", xlabel="Simulated \u03C1")

    # ax.set_title("MPRR Simulated (scaled_rho_sim) vs Estimated Rho (mean)")

    # ax.figure.savefig("MPRR_results.png", dpi=500)

    ax.savefig("rhometa_subsample_results_mixed_sam.png", dpi=500)

    #### deviation plot

    collected_results_sweep_1_df["deviation"] = (collected_results_sweep_1_df["mean"] -
                                                 collected_results_sweep_1_df["scaled_rho_sim"]) / collected_results_sweep_1_df[
                                                    "scaled_rho_sim"]

    ax2 = sns.catplot(data=collected_results_sweep_1_df, x="scaled_rho_sim", y="deviation", hue="Genomes",
                     col="fold_coverage_scaler", col_wrap=2, sharex=True, sharey=True, kind="box", palette="Accent")

    ax2.set(ylim=(-1,1), xlabel="Simulated \u03C1", ylabel="Deviation")

    ax2.map(plt.axhline, y=0, ls='--', color='g', linewidth=2)

    ax2.savefig("rhometa_subsample_deviation_mixed_sam.png", dpi=500)
