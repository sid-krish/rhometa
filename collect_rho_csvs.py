import glob
import pandas as pd

path = './Rho_Est_Output/rho_estimate/'
csv_files = glob.glob(path + "*estimate.csv")

seed_list = [str(file).split('/')[-1].split('_')[0] for file in csv_files]
id_list = [str(file).split('/')[-1].split('_')[1] for file in csv_files]
df_list = [pd.read_csv(file) for file in csv_files]

combined_df = pd.concat(df_list, ignore_index=True)
combined_df["seed"] = seed_list
combined_df["identifier"] = id_list

combined_df = combined_df.reindex(columns=['identifier', 'seed', 'rho','log_likelihood_sum'])
combined_df.to_csv("combined.csv",index=False)