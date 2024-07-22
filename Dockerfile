FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y libgsl-dev && apt-get autoremove && apt-get clean

RUN conda update conda -y

RUN conda install -c bioconda -c conda-forge -y \
  'python=3.12' \
  'numpy=1.26' \
  'pandas=2.2' \
  'scipy=1.14' \
  'seaborn=0.13' \
  'numba=0.60' \
  'tqdm=4.66' \
  'bcftools=1.18' \
  'bwa-mem2=2.2.1' \
  'pysam=0.22' \
  'fastp=0.23.4' \
  'freebayes=1.3.6' \
  'samtools=1.18' \
  'openssl' \ 
  'ncurses' && conda clean -a -y
