import pandas as pd

df = pd.read_csv("10-200_lk_table.csv", nrows=100000, index_col="00 01 10 11")