FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get dist-upgrade --yes
RUN apt-get install -y procps

RUN conda update conda --yes
RUN conda install -c defaults -c bioconda -c conda-forge python pandas numpy seaborn scipy future biopython numba lofreq samtools pysam --yes
RUN conda clean --all --yes