################################################################
## This makefile downloads data for Sequanix demo (#ref)
##
## Author: Claire Rioualen
## Date: 2017-05-23
export LC_ALL=C
export LANG=C

.PHONY: \
	init\
	create_dir\
	download_gene-regulation\
	download_genome_data\
	download_raw_data

#make target FOO=bar

create_dir:
	@echo "Creating ANALYSIS_DIR directory"
	mkdir -p $(ANALYSIS_DIR)


### Download Gene-regulation v4.0

download_gene-regulation:
	cd $(ANALYSIS_DIR) && \
	wget --no-clobber https://github.com/rioualen/gene-regulation/archive/4.0.tar.gz && \
	tar xvzf 4.0.tar.gz && \
	ln -s gene-regulation-4.0 gene-regulation


### Download genome & annotations 

download_genome_data:
	cd $(ANALYSIS_DIR) && \
	wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.30.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome && \
	wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gff3.gz -P ${ANALYSIS_DIR}/genome && \
	wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gtf.gz -P ${ANALYSIS_DIR}/genome && \
	gunzip ${ANALYSIS_DIR}/genome/*.gz


### Download ChIP-seq data 

download_raw_data:
	cd $(ANALYSIS_DIR) && \
	wget --no-clobber ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/005/SRR1176905/SRR1176905.fastq.gz -P ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334674 && \
	wget --no-clobber ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/007/SRR1176907/SRR1176907.fastq.gz -P ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334676 && \
	wget --no-clobber ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/008/SRR1176908/SRR1176908.fastq.gz -P ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334679 && \
	wget --no-clobber ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/000/SRR1176910/SRR1176910.fastq.gz -P ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334677 && \
	gunzip -c ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334674/SRR1176905.fastq.gz > ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334674/GSM1334674.fastq; rm -f ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334674/SRR1176905.fastq.gz && \
	gunzip -c ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334676/SRR1176907.fastq.gz > ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334676/GSM1334676.fastq; rm -f ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334676/SRR1176907.fastq.gz && \
	gunzip -c ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334679/SRR1176908.fastq.gz > ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334679/GSM1334679.fastq; rm -f ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334679/SRR1176908.fastq.gz && \
	gunzip -c ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334677/SRR1176910.fastq.gz > ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334677/GSM1334677.fastq; rm -f ${ANALYSIS_DIR}/ChIP-seq_GSE55357/fastq/GSM1334677/SRR1176910.fastq.gz

all: \
	init\
	create_dir\
	download_gene-regulation\
	download_genome_data\
	download_raw_data

