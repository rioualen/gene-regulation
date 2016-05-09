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

RUN mkdir -p /usr/bin
WORKDIR /usr/bin

RUN git clone https://github.com/rioualen/gene-regulation.git
WORKDIR gene-regulation

RUN make -f scripts/makefiles/install_tools_and_libs_docker.mk all

RUN pip3 install -U pandas

RUN echo export SHELL=/bin/bash >> ~/.bashrc

WORKDIR /usr/bin/gene-regulation

MAINTAINER Claire Rioualen <claire.rioualen@inserm.fr>

