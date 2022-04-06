FROM continuumio/miniconda3:latest

RUN apt-get install -y procps

RUN conda install -c defaults -c bioconda -c conda-forge -y \
  python>=3.9 \
  numpy=1.21 \
  pandas \
  scipy \
  seaborn \
  numba \
  tqdm \
  pysam>=0.17 \
  freebayes>=1.3.5 \
  samtools==1.15 \
  openssl=1.1.1n