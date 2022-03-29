FROM continuumio/miniconda3

RUN apt-get install -y procps

RUN conda install -c defaults -c bioconda -c conda-forge -y \
    numpy \
    pandas \
    scipy \
    pysam \
    seaborn \
    numba \
    samtools=1.12 \
    freebayes>=1.3.5