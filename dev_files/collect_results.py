import glob

import pandas as pd
import seaborn as sns


def theta_results(theta_results_path):
    csv_files = glob.glob(theta_results_path + "*.csv")

    sample_list = [str(file).split('/')[-1].split('_')[0] for file in csv_files]
    id_list = [str(file).split('/')[-1].split('_')[1] for file in csv_files]

    df_list = [pd.read_csv(file, skiprows=2, header=None).T.drop(index=0) for file in csv_files]
    combined_df = pd.concat(df_list, ignore_index=True)

    cols = ["mean_depth", "tps_mean_depth", "median_depth", "tps_median_depth"]
    combined_df.columns = cols
    combined_df["sample"] = sample_list
    combined_df["identifier"] = id_list

    combined_df = combined_df.reindex(
        columns=["sample", "identifier", "mean_depth", "tps_mean_depth", "median_depth", "tps_median_depth"])
    combined_df.to_csv("collected_theta_results.csv", index=False)

    sns.set_palette("Set1")
    fig = sns.boxplot(combined_df, x="tps_mean_depth")
    fig.figure.savefig("theta_results_plot.png", dpi=500)

    return combined_df


def rho_results(rho_results_path):
    csv_files = glob.glob(rho_results_path + "*estimate.csv")

    sample_list = [str(file).split('/')[-1].split('_')[0] for file in csv_files]
    seed_list = [str(file).split('/')[-1].split('_')[1] for file in csv_files]
    id_list = [str(file).split('/')[-1].split('_')[2] for file in csv_files]

    df_list = [pd.read_csv(file) for file in csv_files]

    combined_df = pd.concat(df_list, ignore_index=True)

    combined_df["sample"] = sample_list
    combined_df["seed"] = seed_list
    combined_df["identifier"] = id_list

    combined_df = combined_df.reindex(columns=['sample', 'identifier', 'seed', 'rho', 'log_likelihood_sum'])
    combined_df.to_csv("collected_rho_results.csv", index=False)

    sns.set_palette("Set2")
    fig = sns.boxplot(combined_df, x="rho")
    fig.figure.savefig("rho_results_plot.png", dpi=500)

    return combined_df


if __name__ == '__main__':
    theta_results_path = './Theta_Est_Output/theta_estimate/'
    rho_results_path = './Rho_Est_Output/rho_estimate/'

    theta_combined_df = theta_results(theta_results_path)
    rho_combined_df = rho_results(rho_results_path)

    merged_df = pd.merge(rho_combined_df, theta_combined_df, on='identifier').drop(columns=["sample_y"])
    merged_df.to_csv("collected_merged_results.csv", index=False)