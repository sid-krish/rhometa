FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get upgrade -y

RUN apt-get install -y procps && apt-get clean -y

RUN conda update conda -y

RUN conda install -c bioconda -c conda-forge -y \
  'python=3.10' \
  'numpy=1.23' \
  'pandas=1.5' \
  'scipy=1.9' \
  'seaborn=0.12' \
  'numba=0.56' \
  'tqdm=4.64' \
  'bcftools=1.16' \
  'bwa-mem2=2.2.1' \
  'pysam=0.20' \
  'fastp=0.23.2' \
  'freebayes=1.3.6' \
  'samtools=1.16' \
  'openssl>=1.1.1t'\ 
  'ncurses>=6.4' && conda clean -a -y
