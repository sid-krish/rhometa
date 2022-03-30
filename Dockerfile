FROM continuumio/miniconda3:latest

# update the base system
RUN apt-get update && \
    apt-get dist-upgrade --yes && \
    apt-get install -y procps && \
    apt-get clean

# just install non-python deps with conda
RUN conda update conda --yes && \
    conda install --yes  -c conda-forge -c bioconda \
        bcftools \
        bwa \
        "freebayes>=1.3.6" \
        "msprime>=1.1.0" \
        python \
        "samtools>=1.15" && \
    conda clean --all

# install all python modules using pip
RUN pip install --no-cache-dir \
        biopython \
        future \
        scipy \
        seaborn \
        tqdm \
        pandas \
        tables \
        pysam \
        numba \
        numpy
