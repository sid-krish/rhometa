FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get dist-upgrade --yes
RUN apt-get install samtools procps --yes

RUN conda update conda --yes
RUN conda install -c main -c bioconda pandas seaborn scikit-learn python=3.6 scipy future biopython numba plotly fastsimbac seq-gen bwa lofreq pysam cyvcf2 --yes
RUN pip install pyarrow
RUN conda clean --all --yes