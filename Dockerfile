FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get dist-upgrade --yes
RUN apt-get install -y procps

RUN conda update conda --yes
RUN conda install -c defaults -c bioconda -c conda-forge python pandas seaborn scikit-learn scipy future biopython numba fastsimbac seq-gen bwa lofreq cyvcf2 samtools pysam --yes
RUN conda clean --all --yes