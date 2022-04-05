FROM continuumio/miniconda3:latest

RUN apt-get install -y procps

RUN conda install -c defaults -c bioconda -c conda-forge -y \
  python \
  numpy=1.21 \
  pandas \
  scipy \
  pysam \
  seaborn \
  numba \
  tqdm \
  freebayes=1.3.5 \
  samtools==1.12 \
  bcftools==1.12