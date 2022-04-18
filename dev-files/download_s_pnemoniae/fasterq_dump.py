from sys import argv
import subprocess

accessions_txt = argv[1]

with open(accessions_txt, 'r') as f:
    accessions = f.readlines()
    accessions = [i.strip() for i in accessions]
    for i in accessions:
        print(i)
        subprocess.call(f"fasterq-dump --split-3 --threads 4 {i}", shell=True)
