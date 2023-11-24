FROM continuumio/miniconda3

LABEL "docker.omicstp"="MasterBBC"
LABEL version="1.0"

# mettre a jour conda et install python 3.8 et les outils par conda
RUN conda update -n base -c defaults conda && \
	conda install -c bioconda -y python=3.8 bowtie2 subread fastqc && \
	conda clean --all --yes

# installation des outils en pip 
# (leur installation en conda ne marche pas)
RUN pip install multiqc cutadapt && \
	pip cache purge

# samtools
RUN apt-get update && \
	apt-get -y install wget unzip bzip2 curl gcc make && \
    apt-get -y install libz-dev libncurses5-dev && \
    apt-get -y install libbz2-dev liblzma-dev && \
    wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
	tar xvjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --prefix /usr/bin/samtools-1.18 && \
    make && \
    make install && \
	apt-get -qq -y autoremove && \
  	apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ENV PATH=$PATH:/usr/bin/samtools-1.18/bin




    
 
