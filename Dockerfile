#############################################################
# Dockerfile to build a running env container for starscope #
#############################################################

## From conda docker
FROM continuumio/miniconda3

## Maintainer
MAINTAINER oben <obennoname@gmail.com>

## Setup conda env
USER root
## install tini for init process
RUN apt update && apt install -y git tini
## use mamba instead of conda command
RUN conda install mamba -n base -c conda-forge

## Setup workdir
WORKDIR /app

## download app source
##RUN git clone https://github.com/obenno/scRNA-seq ./scRNA-seq
##WORKDIR /app/scRNA-seq
##COPY scRNAseq_env.yml .
## create conda env with app env file
##RUN mamba env create -f scRNAseq_env.yml

## All packages will be installed to base env
## and could be invoked directoryly
RUN mamba install -c conda-forge -c bioconda plotly samtools==1.15 star==2.7.11a cutadapt fastx_toolkit fastqc r-flexdashboard r-rmarkdown r-scales r-tidyverse r-plotly r-kableextra r-jsonlite r-dt r-future r-seurat r-seuratdisk r-rhpcblasctl bedtools bioconductor-limma gawk pigz gzip coreutils jq qualimap

## install SeuratDisk
##RUN mamba install -c bioconda -c conda-forge r-hdf5r r-cli r-crayon r-matrix r-r6 r-rlang r-withr r-stringi
##RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes", repos = "https://cloud.r-project.org/", lib = "/opt/conda/lib/R/library")}; remotes::install_github("mojaveazure/seurat-disk", lib = "/opt/conda/lib/R/library")'
##
#### install RhpcBLASctl to control BLAS threads
##RUN mamba install r-rhpcblasctl

## install trust4 for VDJ analysis
RUN mamba install -c conda-forge -c bioconda trust4

## copy entrypoint.sh
##COPY entrypoint.sh .

## Follow Dockstore's guide
## switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -m -r -g ubuntu -u 1000 ubuntu
RUN chown -R ubuntu: /app
USER ubuntu

## Setup bashrc file for the ubuntu user
##RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
##RUN echo "conda activate starscope_env" >> ~/.bashrc

##SHELL ["/bin/bash", "--login", "-c"]
## setup default cmd (deprecated)
## no entrypoint and cmd will be set
##RUN chmod +x entrypoint.sh
##ENTRYPOINT ["./entrypoint.sh"]
