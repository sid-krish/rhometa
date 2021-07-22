FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get dist-upgrade --yes
RUN apt-get install -y procps

RUN conda update conda --yes
RUN conda install -c bioconda -c conda-forge pandas <= 1.2 seaborn scikit-learn scipy future biopython numba plotly fastsimbac seq-gen bwa lofreq cyvcf2 samtools pysam pyarrow --yes
RUN conda clean --all --yes