FROM continuumio/miniconda3

LABEL "docker.omicstp"="MasterBBC"
LABEL version="1.0"

# mettre a jour conda et install python 3.8 et les outils par conda
RUN conda update -n base -c defaults conda && \
	conda install -c bioconda -y python=3.8 bowtie2 subread samtools fastqc && \
	conda clean --all --yes

# installation des outils en pip 
# (leur installation en conda ne marche pas)
RUN pip install multiqc cutadapt && \
	pip cache purge

    
