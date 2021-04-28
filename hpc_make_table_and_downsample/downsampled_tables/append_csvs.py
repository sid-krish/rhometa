import glob

import pandas as pd

# File must be in and executed from working directory

extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

all_filenames.sort(key=lambda x: int(x[:-4].split('_')[2]))

all_filenames = all_filenames[0:101] # for testing just a few files

# combine all files in the list
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])

# export to csv
combined_csv.to_csv("10-100_lk_table.csv", index=False)
