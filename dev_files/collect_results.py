import glob

import numpy as np
import pandas as pd
import seaborn as sns


def cov_results(cov_results_path):
    cov_files = glob.glob(cov_results_path + "*.cov")

    reference_list = [str(file).split('/')[-1].split('_')[-2]
                      for file in cov_files]
    id_list = [str(file).split('/')[-1].split('_')[0] for file in cov_files]

    df_list = [pd.read_table(file, usecols=["coverage", "endpos"])
               for file in cov_files]

    weighted_mean_cov_list = []
    for i in df_list:
        coverage = i["coverage"].to_numpy(dtype="float32")
        genome_lens = i["endpos"].to_numpy(dtype="int32")

        weight_covs = genome_lens * coverage
        sum_weight_covs = np.sum(weight_covs)
        genome_lens_sum = np.sum(genome_lens)

        weighted_mean_cov = sum_weight_covs / genome_lens_sum
        weighted_mean_cov_list.append(weighted_mean_cov)

    combined_df = pd.DataFrame()
    combined_df["reference"] = reference_list
    combined_df["identifier"] = id_list
    combined_df["weighted_mean_cov"] = weighted_mean_cov_list

    combined_df = combined_df.reindex(
        columns=["reference", "identifier", "weighted_mean_cov"])
    combined_df.to_csv("collected_coverage_results.csv", index=False)

    sns.set_palette("Set1")
    fig = sns.boxplot(combined_df, x="weighted_mean_cov")
    fig.set_xlim(1, 100)
    fig.figure.savefig("coverage_results_plot.png", dpi=500)
    fig.figure.clf()

    return combined_df


def theta_results(theta_results_path):
    csv_files = glob.glob(theta_results_path + "*.csv")

    reference_list = [str(file).split('/')[-1].split('_')[-5]
                      for file in csv_files]
    id_list = [str(file).split('/')[-1].split('_')[0] for file in csv_files]

    df_list = [pd.read_csv(file, skiprows=2, header=None).T.drop(
        index=0) for file in csv_files]
    combined_df = pd.concat(df_list, ignore_index=True)

    cols = ["mean_depth", "tps_mean_depth", "median_depth", "tps_median_depth"]
    combined_df.columns = cols
    combined_df["reference"] = reference_list
    combined_df["identifier"] = id_list

    combined_df = combined_df.reindex(
        columns=["reference", "identifier", "mean_depth", "tps_mean_depth", "median_depth", "tps_median_depth"])
    combined_df.to_csv("collected_theta_results.csv", index=False)

    sns.set_palette("Set2")
    fig = sns.boxplot(combined_df, x="tps_mean_depth")
    fig.figure.savefig("theta_results_plot.png", dpi=500)
    fig.figure.clf()

    return combined_df


def rho_results(rho_results_path):
    csv_files = glob.glob(rho_results_path + "*estimate.csv")

    reference_list = [str(file).split('/')[-1].split('_')[-4]
                      for file in csv_files]
    seed_list = [str(file).split('/')[-1].split('_')[0] for file in csv_files]
    id_list = [str(file).split('/')[-1].split('_')[1] for file in csv_files]

    df_list = [pd.read_csv(file) for file in csv_files]

    combined_df = pd.concat(df_list, ignore_index=True)

    combined_df["reference"] = reference_list
    combined_df["seed"] = seed_list
    combined_df["identifier"] = id_list

    combined_df = combined_df.reindex(
        columns=['reference', 'identifier', 'seed', 'rho', 'log_likelihood_sum'])
    combined_df.to_csv("collected_rho_results.csv", index=False)

    sns.set_palette("Set3")
    fig = sns.boxplot(combined_df, x="rho")
    fig.set_xlim(1, 100)
    fig.figure.savefig("rho_results_plot.png", dpi=500)
    fig.figure.clf()

    return combined_df


if __name__ == '__main__':
    cov_results_path = ''
    theta_results_path = ''
    rho_results_path = ''

    cov_combined_df = cov_results(cov_results_path)
    theta_combined_df = theta_results(theta_results_path)
    rho_combined_df = rho_results(rho_results_path)

    cov_theta_merged_df = pd.merge(cov_combined_df, theta_combined_df,
                                   on=['reference', 'identifier'])

    cov_theta_rho_merged_df = pd.merge(cov_theta_merged_df, rho_combined_df,
                                       on=['reference', 'identifier'])

    cov_theta_rho_merged_df = cov_theta_rho_merged_df.sort_values(
        by=['reference', 'identifier'])

    cov_theta_rho_merged_df.to_csv("collected_merged_results.csv", index=False)
