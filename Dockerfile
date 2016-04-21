FROM ubuntu:14.04

FROM adreeve/python-numpy:latest

RUN apt-get update \
    && apt-get --yes install git \
    wget \
    nano \
    python-software-properties\
    build-essential \
    python3-software-properties \
    software-properties-common

#RUN mkdir -p /data/GSE20870/GSM521934 /data/GSE20870/GSM521935
#WORKDIR /data/GSE20870/GSM521934
#RUN wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra
#WORKDIR /data/GSE20870/GSM521935
#RUN wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra

RUN mkdir -p /usr/bin
WORKDIR /usr/bin

RUN git clone https://github.com/rioualen/gene-regulation.git
WORKDIR gene-regulation

## todo replace makefile with docker images...
RUN make -f scripts/makefiles/install_tools_and_libs_docker.mk all



## todo upgrade pandas
RUN pip3 install -U pandas

RUN echo export SHELL=/bin/bash >> ~/.bashrc

WORKDIR /usr/bin/gene-regulation

## todo as bin bash in bashrc

## todo add entry point

MAINTAINER Claire Rioualen <claire.rioualen@inserm.fr>

## usage
## cd gene-regulation
## docker build -t gene-regulation/test:latest .
## docker run -v /data/results/GSE20870_/:/data/results/GSE20870/ -it gene-regulation/test:latest . /bin/bash
## docker run -v /data/results/GSE20870_/:/data/results/GSE20870/ -v /data/GSE20870/:/data/GSE20870/ -it rioualen/gene-regulation:0.3 /bin/bash
## docker push rioualen/gene-regulation:0.3

