#############################################################
# Dockerfile to build a running env container for starscope #
#############################################################

## From official miniforge image
## https://github.com/conda-forge/miniforge-images
FROM quay.io/condaforge/miniforge-pypy3

## Maintainer
MAINTAINER oben <obennoname@gmail.com>

## Setup conda env
USER root
## install tini for init process
RUN apt update && apt install -y git tini

## Setup workdir
WORKDIR /app

## create conda env with app env file
COPY scRNAseq_env.yml .
RUN mamba env create -f scRNAseq_env.yml

## copy entrypoint.sh
COPY entrypoint.sh .


## Setup bashrc file for the ubuntu user
##RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
##RUN echo "conda activate starscope_env" >> ~/.bashrc
##SHELL ["/bin/bash", "--login", "-c"]
## setup default cmd (deprecated)
## no entrypoint and cmd will be set
RUN chmod +x entrypoint.sh

## Follow Dockstore's guide
## switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -m -r -g ubuntu -u 1000 ubuntu
RUN chown -R ubuntu: /app
USER ubuntu

ENTRYPOINT ["./entrypoint.sh"]
