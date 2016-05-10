FROM ubuntu:14.04

FROM rioualen/ubuntu-essentials:1.0
FROM rioualen/python-essentials:1.0
FROM rioualen/r-essentials:1.0

FROM rioualen/fastqc:0.11.5
FROM rioualen/sickle:1.33
FROM rioualen/bedtools:2.24.0
FROM rioualen/samtools:1.3
FROM rioualen/bowtie:1.1.1




FROM rioualen/sratoolkit:2.5.2
FROM rioualen/bwa:latest




FROM rioualen/macs2:latest
FROM rioualen/macs:1.4.2
FROM rioualen/homer:latest
FROM rioualen/spp:1.11
FROM rioualen/bowtie2:2.2.6


ENV WORKSPACE=/root/workspace
WORKDIR ${WORKSPACE}

RUN wget https://github.com/rioualen/gene-regulation/archive/1.0.tar.gz
RUN tar zvxf 1.0.tar.gz
WORKDIR gene-regulation-1.0

#RUN make -f scripts/makefiles/install_tools_and_libs_docker.mk all


RUN echo export SHELL=/bin/bash >> /root/.bashrc
RUN echo export PATH=$PATH:/root/bin >> /root/.bashrc

ENTRYPOINT ["/bin/bash"]


MAINTAINER Claire Rioualen <claire.rioualen@inserm.fr>

