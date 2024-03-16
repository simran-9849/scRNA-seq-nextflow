#############################################################
# Dockerfile to build a running env container for starscope #
#############################################################

## From official miniforge image
## https://github.com/conda-forge/miniforge-images
FROM quay.io/condaforge/miniforge3

## Maintainer
MAINTAINER Zhixia Xiao <obennoname@gmail.com>

## Setup conda env
USER root
## install tini for init process
## tini was already used by miniforge-pypy3
##RUN apt update && apt install -y git tini
## register Tini as a child subreaper
ENV TINI_SUBREAPER=true

## Setup workdir
WORKDIR /app

## Update the base env with app env file
## All the packages will be installed to base env and no need to activate
## specific package in entrypoint.sh
COPY scRNAseq_env.yml .
RUN mamba env update -n base -f scRNAseq_env.yml && rm scRNAseq_env.yml

## copy entrypoint.sh
##COPY entrypoint.sh .

##RUN chmod +x entrypoint.sh

## Follow Dockstore's guide
## switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -m -r -g ubuntu -u 1000 ubuntu
RUN chown -R ubuntu: /app
USER ubuntu

##ENTRYPOINT ["/app/entrypoint.sh"]
