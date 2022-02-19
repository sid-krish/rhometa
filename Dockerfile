FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get dist-upgrade --yes
RUN apt-get install -y procps

RUN conda update conda --yes
RUN conda install -c defaults -c bioconda -c conda-forge python numpy openblas pandas seaborn scipy future biopython numba fastsimbac seq-gen bwa lofreq pysam samtools=1.12 --yes
RUN conda install -c conda-forge msprime>=1.1.0 --yes
RUN conda clean --all --yes