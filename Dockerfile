#############################################################
# Dockerfile to build a running env container for starscope #
#############################################################

## From conda docker
FROM condaforge/miniforge-pypy3

## Maintainer
MAINTAINER oben <obennoname@gmail.com>

## Setup conda env
USER root
## install tini for init process
## tini was already used by miniforge-pypy3
##RUN apt update && apt install -y git tini
## use mamba instead of conda command
##RUN conda install mamba -n base -c conda-forge

## Setup workdir
WORKDIR /app

COPY scRNAseq_env.yml .
## https://stackoverflow.com/questions/42352841/how-to-update-an-existing-conda-environment-with-a-yml-file
## update conda base env with app env file
##RUN mamba env update -n base -f scRNAseq_env.yml --prune
RUN mamba env create -f scRNAseq_env.yml && mamba clean --force-pkgs-dirs --all --yes

## Setup bashrc file for the all the user
RUN sed -i 's/conda activate base/conda activate starscope_env/' ~/.bashrc && sed -i 's/conda activate base/conda activate starscope_env/' /etc/skel/.bashrc

COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

## Follow Dockstore's guide
## switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -m -r -g ubuntu -u 1000 ubuntu
RUN chown -R ubuntu: /app
USER ubuntu

##SHELL ["/bin/bash", "--login", "-c"]
## setup default cmd (deprecated)
## no entrypoint and cmd will be set

## Add ENTRYPOINT to ensure user does not exist could activate the default conda environment
ENTRYPOINT ["/app/entrypoint.sh"]
