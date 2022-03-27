FROM continuumio/miniconda3

RUN apt-get update && apt-get dist-upgrade --yes
RUN apt-get install -y procps

RUN conda update conda --yes
RUN conda install -c defaults -c bioconda -c conda-forge \
    python \
    numpy>=1.22 \
    pandas>=1.4.0 \
    biopython \
    scipy \
    future \
    bwa \
    samtools>=1.15 \
    pysam>=0.18.0 \
    seaborn \
    numba \
    msprime>=1.1.0 \
    tqdm \
    bcftools \
    freebayes>=1.3.6 \
    libopenblas --yes
RUN conda clean --all --yes