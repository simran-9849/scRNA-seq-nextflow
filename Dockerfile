#############################################################
# Dockerfile to build a sample tool container for shinysyn
#############################################################

## From conda docker
FROM continuumio/miniconda3

## Maintainer
MAINTAINER oben <obennoname@gmail.com>

## Setup conda env
USER root
RUN apt update && apt install git
## use mamba instead of conda command
RUN conda install mamba -n base -c conda-forge

## Setup workdir
WORKDIR /app

## download app source
##RUN git clone https://github.com/obenno/scRNA-seq ./scRNA-seq
WORKDIR /app/scRNA-seq
COPY scRNAseq_env.yml .
## create conda env with app env file
RUN mamba env create -f scRNAseq_env.yml

## copy entrypoint.sh
COPY entrypoint.sh .

## Follow Dockstore's guide
## switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -m -r -g ubuntu -u 1000 ubuntu
RUN chown -R ubuntu: /app/scRNA-seq
USER ubuntu

# Setup bashrc file for the ubuntu user
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate scRNAseq_env" >> ~/.bashrc

##SHELL ["/bin/bash", "--login", "-c"]
## setup default cmd
RUN chmod +x entrypoint.sh
ENTRYPOINT ["./entrypoint.sh"]
