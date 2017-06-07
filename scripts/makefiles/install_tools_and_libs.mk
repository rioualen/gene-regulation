################################################################
## This makefile installs bioinformatics tools for the analysis of NGS
## data about gene regulation, together with the required libraries.
##
## It is conceived for Linux Ubuntu systems mlanaged by apt-get. The
## targets can probably be adapted to other systems as well.
##
## By default, tools are installed in a user-custom directory
## $$HOME/bin. This can be changed by modifying the USER_BIN
## directory.
##
## Author: Claire Rioualen
## Date: 2016

export LC_ALL=C
export LANG=C

.PHONY: \
	create_bin \
	init_BASHRC \
	init \
	add_pub_key \
	add_repos \
	add_packages \
	R_installation \
	R_lib \
	python \
	java \
	samtools \
	bedtools \
	sratoolkit \
	fastqc \
	sickle \
	bowtie \
	bowtie2 \
	bwa \
	spp \
	macs1 \
	macs2 \
	homer \
	swembl \
	igv \
	igv_tools \
	desktop_and_x2go \
	Rstudio 


################################################################
## List targets
MAKEFILE=scripts/makefiles/install_tools_and_libs.mk
usage:
	echo "usage: make [-OPT='options'] target"
	echo "implemented targets"
	perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

## Tool & versions


CRAN_PACK_LIST='XML', 'bPeaks', 'caTools', 'VennDiagram', 'devtools'
BIOC_PACK_LIST='affy', 'biomaRt', 'Rsamtools', 'genefilter', 'GenomicFeatures', 'edgeR', 'DESeq2'
PUB_KEY="51716619E084DAB9"
#E084DAB9 51716619E084DAB9 F7B8CEA6056E8E56 06F90DE5381BA480

BEDOPS_VER=2.4.19
BEDTOOLS_VER=2.24.0
BOWTIE1_VER=1.1.1
BOWTIE2_VER=2.2.6
CUFFLINKS_VER=2.2.1
DEEPTOOLS_VER=2.5.0.1
FASTQC_VER=0.11.5
IGV_VER=2.3.59
IGVTOOLS_VER=2.3.57
JAVA_VER=8
MACS1_VER=1.4.2
NUMPY_VER=1.9.2
RPY2_VER=2.5.6
RSTUDIO_VER=0.99.903
SAMTOOLS_VER=1.3
SRATOOLKIT_VER=2.5.2
SCIPY_VER=0.16.0
SICKLE_VER=1.33
SPP_VER=1.11
STAR_VER=2.5.2b
SUBREAD_VER=1.5.2
SWEMBL_VER=3.3.1

# ================================================================
# General configuration
# ================================================================

# ----------------------------------------------------------------
# bin & bashrc
# ----------------------------------------------------------------

## Directory where all the executables will be installed
BIN_DIR=$(HOME)/bin

## Directory in which the source code will be downloaded and compiled when required
SOURCE_DIR=$(HOME)/app_sources

BASHRC=$(HOME)/.bashrc

## Create the bin & source directories
create_bin:
	@echo "Creating BIN_DIR and SOURCE_DIR directories"
	mkdir -p $(BIN_DIR)
	mkdir -p $(SOURCE_DIR)
	@echo "	BIN_DIR	${BIN_DIR}"
	@echo "	SOURCE_DIR	${SOURCE_DIR}"

## Initialize a bashrc file, which will allow any user to load the
## path and configuration at once.
update_bashrc:
	@echo 'alias ls="ls --color"' >> $(BASHRC)
	@echo 'alias rm="rm -i"' >> $(BASHRC)
	@echo "" >> $(BASHRC)
	@echo "export LC_ALL=C" >> $(BASHRC)
	@echo "export LANG=C" >> $(BASHRC)
	@echo 'export force_color_prompt=yes' >> $(BASHRC)
	@echo 'export LC_COLLATE=C' >> $(BASHRC)
	@echo "" >> $(BASHRC)
	@echo 'export PATH='$(PATH):$(BIN_DIR) >> $(BASHRC)

init: create_bin update_bashrc

# ----------------------------------------------------------------
# Basic libraries: Ubuntu, R, python, java
# ----------------------------------------------------------------

## ???
add_pub_key:
	for i in '$(PUB_KEY)'; do echo "PUB_KEY: $$i"; sudo apt-key adv --recv-keys --keyserver keyserver.ubuntu.com $$i; done

add_repos: 
	sudo apt-get install --yes python-software-properties
	sudo add-apt-repository ppa:x2go/stable --yes
	sudo apt-add-repository ppa:ubuntu-mate-dev/ppa --yes
	sudo apt-add-repository ppa:ubuntu-mate-dev/trusty-mate --yes
	sudo add-apt-repository ppa:ravefinity-project/ppa --yes
	sudo apt-get update --yes

add_packages:
	sudo apt-get update
	sudo apt-get -y install ssh rsync git graphviz gedit-plugins  apt-transport-https
	sudo apt-get -y install zlibc zlib1g-dev									# Required by sickle, bamtools, samtools...
	sudo apt-get -y install build-essential										# Includes gcc compiler
	sudo apt-get -y install libncurses5-dev libncursesw5-dev					# Required at least by samtools
	sudo apt-get -y install libboost-dev										# might be required for spp, to be checked
	sudo apt-get -y install gdebi												# required by rstudio install

R_installation:

### !!! sudo bash -c "echo 'deb http://cran.univ-lyon1.fr/bin/linux/ubuntu trusty/' >> /etc/apt/sources.list"

	sudo echo "deb http://cran.univ-lyon1.fr/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
	sudo apt-get update
	sudo apt-get --force-yes install r-base r-base-dev libcurl4-openssl-dev libxml2-dev
	echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" >> ~/.Rprofile



R_lib: 
	sudo Rscript -e "pack.list <- c($(CRAN_PACK_LIST));\
	pack <- pack.list[!(pack.list %in% installed.packages()[,'Package'])];\
	if(length(pack)) install.packages(pack);\
	source('http://bioconductor.org/biocLite.R');\
	pack.list <- c($(BIOC_PACK_LIST));\
	pack <- pack.list[!(pack.list %in% installed.packages()[,'Package'])];\
	if(length(pack)) biocLite(pack);\
	library(devtools);\
	install_github('PF2-pasteur-fr/SARTools', build_vignettes=TRUE)"

Rstudio: 
	cd $(SOURCE_DIR) && \
	sudo apt-get install -y libjpeg62
	wget --no-clobber https://download1.rstudio.org/rstudio-$(RSTUDIO_VER)-amd64.deb
	yes | sudo gdebi rstudio-$(RSTUDIO_VER)-amd64.deb
	rm -f rstudio-$(RSTUDIO_VER)-amd64.deb

python: 
	sudo apt-get -y install python-pip python-dev
	sudo apt-get -y install python3-pip python3.4-dev
	sudo pip3 install "rpy2<$(RPY2_VER)"
	sudo pip3 install numpy
	sudo pip install numpy
	sudo pip3 install snakemake docutils pandas
	sudo pip3 install pyyaml
	sudo pip3 install psutil
	sudo pip install pyBigWig
	sudo pip install pysam==0.8.4


ubuntu_packages: add_repos add_packages

R: R_installation R_lib Rstudio


# ================================================================
# NGS
# ================================================================

# ----------------------------------------------------------------
# General tools
# ----------------------------------------------------------------

samtools:
	cd $(SOURCE_DIR) && \
	wget --no-clobber http://sourceforge.net/projects/samtools/files/samtools/$(SAMTOOLS_VER)/samtools-$(SAMTOOLS_VER).tar.bz2 && \
	bunzip2 -f samtools-$(SAMTOOLS_VER).tar.bz2 && \
	tar xvf samtools-$(SAMTOOLS_VER).tar && \
	cd samtools-$(SAMTOOLS_VER) && \
	make  && \
	sudo make install

bedtools:
	cd $(SOURCE_DIR) && \
	wget --no-clobber https://github.com/arq5x/bedtools2/releases/download/v$(BEDTOOLS_VER)/bedtools-$(BEDTOOLS_VER).tar.gz && \
	tar xvfz bedtools-$(BEDTOOLS_VER).tar.gz && \
	cd bedtools2 && \
	make && \
	sudo make install

sratoolkit:
	cd $(SOURCE_DIR) && \
	wget -nc http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$(SRATOOLKIT_VER)/sratoolkit.$(SRATOOLKIT_VER)-ubuntu64.tar.gz && \
	tar xzf sratoolkit.$(SRATOOLKIT_VER)-ubuntu64.tar.gz && \
	cp `find sratoolkit.$(SRATOOLKIT_VER)-ubuntu64/bin -maxdepth 1 -executable -type l` $(BIN_DIR)

bedops:
	cd $(SOURCE_DIR) && \
	wget -nc https://github.com/bedops/bedops/releases/download/v$(BEDOPS_VER)/bedops_linux_x86_64-v$(BEDOPS_VER).tar.bz2 && \
	tar jxvf bedops_linux_x86_64-v$(BEDOPS_VER).tar.bz2 && \
	mkdir -p bedops && \
	mv bin bedops && \
	cp bedops/bin/* $(BIN_DIR)

deeptools:
	cd $(SOURCE_DIR) && \
	wget -nc https://github.com/fidelram/deepTools/archive/$(DEEPTOOLS_VER).tar.gz && \
	tar xvzf $(DEEPTOOLS_VER).tar.gz && \
	cd  deepTools-$(DEEPTOOLS_VER) && \
	sudo python setup.py install

## Picard tools todo
#! Java 8 required
#https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar
#create $PICARD env

# ----------------------------------------------------------------
# Quality assessment & trimming
# ----------------------------------------------------------------

fastqc:
	cd $(SOURCE_DIR) && \
	wget --no-clobber http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v$(FASTQC_VER).zip && \
	unzip -o fastqc_v$(FASTQC_VER).zip && \
	chmod +x FastQC/fastqc && \
	ln -s -f $(SOURCE_DIR)/FastQC/fastqc $(BIN_DIR)/fastqc

sickle:
	cd $(SOURCE_DIR) && \
	wget -nc https://github.com/najoshi/sickle/archive/v$(SICKLE_VER).tar.gz && \
	tar xvzf v$(SICKLE_VER).tar.gz && \
	cd sickle-$(SICKLE_VER) && \
	make && \
	cp sickle $(BIN_DIR)

# ----------------------------------------------------------------
# Mapping tools
# ----------------------------------------------------------------

bowtie: 
	cd $(SOURCE_DIR) && \
	wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie/$(BOWTIE1_VER)/bowtie-$(BOWTIE1_VER)-linux-x86_64.zip && \
	unzip bowtie-$(BOWTIE1_VER)-linux-x86_64.zip && \
	cp `find bowtie-$(BOWTIE1_VER)/ -maxdepth 1 -executable -type f` $(BIN_DIR)

bowtie2:
	cd $(SOURCE_DIR) && \
	wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$(BOWTIE2_VER)/bowtie2-$(BOWTIE2_VER)-linux-x86_64.zip && \
	unzip bowtie2-$(BOWTIE2_VER)-linux-x86_64.zip && \
	cp `find bowtie2-$(BOWTIE2_VER)/ -maxdepth 1 -executable -type f` $(BIN_DIR)

bwa:
	sudo apt-get -y install bwa
#	wget -nc https://sourceforge.net/projects/bio-bwa/files/bwa-$(BWA_VER).tar.bz2 && \

subread:
	cd $(SOURCE_DIR) && \
	wget -nc https://sourceforge.net/projects/subread/files/subread-$(SUBREAD_VER)/subread-$(SUBREAD_VER)-source.tar.gz && \
	tar zxvf subread-$(SUBREAD_VER)-source.tar.gz && \
	cd subread-$(SUBREAD_VER)-source/src && \
	make -f Makefile.Linux && \
	cd ../bin && \
	cp `find * -executable -type f` $(BIN_DIR)

#star:
#	cd $(SOURCE_DIR) && \
#	wget -nc https://github.com/alexdobin/STAR/archive/$(STAR_VER).tar.gz && \
#	tar xvzf $(STAR_VER).tar.gz && \
#	cd STAR-$(STAR_VER) && \
#	make STAR && \


tophat2:
	cd $(SOURCE_DIR) && \
	wget --no-clobber https://ccb.jhu.edu/software/tophat/downloads/tophat-$(TOPHAT_VER).Linux_x86_64.tar.gz && \
	tar xvfz tophat-$(TOPHAT_VER).Linux_x86_64.tar.gz && \
	cd tophat-$(TOPHAT_VER).Linux_x86_64 && \
	rm -Rf AUTHORS LICENSE README intervaltree/ sortedcontainers/ && \
	mv ./* $(BIN_DIR) && \
	cd ..; rm -Rf tophat-$(TOPHAT_VER).Linux_x86_64*


# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------

#macs1:
#	cd $(SOURCE_DIR) && \
#	wget --no-clobber https://github.com/downloads/taoliu/MACS/MACS-$(MACS1_VER)-1.tar.gz && \
#	tar xvfz MACS-$(MACS1_VER)-1.tar.gz  && \
#	cd MACS-$(MACS1_VER) && \
#	sudo python setup.py install && \

macs1:
	cd $(SOURCE_DIR) && \
	wget --no-clobber https://pypi.python.org/packages/86/da/1e57f6e130b732160d87d96f2cc1771b9de24ce16522a4f73a8528166b87/MACS-1.4.3.tar.gz && \
	tar xvzf MACS-1.4.3.tar.gz && \
	cd MACS-1.4.3 && \
	sudo python setup.py install


macs2:
	sudo pip install MACS2

#spp:
#	cd $(SOURCE_DIR) && \
#	wget -nc http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_$(SPP_VER).tar.gz && \
#	sudo R CMD INSTALL spp_$(SPP_VER).tar.gz

## Currently not working -> see w/ S.Wilder, new version to come?
#swembl:
#	cd $(SOURCE_DIR) && \
#	wget "http://www.ebi.ac.uk/~swilder/SWEMBL/SWEMBL.$(SWEMBL_VER).tar.bz2" && \
#	bunzip2 -f SWEMBL.$(SWEMBL_VER).tar.bz2 && \
#	tar xvf SWEMBL.$(SWEMBL_VER).tar && \
#	rm SWEMBL.$(SWEMBL_VER).tar && \
#	chown -R ubuntu-user SWEMBL.$(SWEMBL_VER) && \
#	cd SWEMBL.$(SWEMBL_VER) && \
#	make
##gcc main.c IO.c calc.c stack.c summit.c refcalc.c wiggle.c overlap.c -o SWEMBL -lz -lm

homer:
	mkdir $(SOURCE_DIR)/homer && \
	cd $(SOURCE_DIR)/homer && \
	wget -nc "http://homer.salk.edu/homer/configureHomer.pl" && \
	perl configureHomer.pl -install homer && \
	cp `find $(SOURCE_DIR)/homer/bin -maxdepth 1 -executable -type f` $(BIN_DIR)

## TODO
#source("https://bioconductor.org/biocLite.R")
#biocLite("mosaics")

##install.packages("devtools")
#library(devtools)
#install_github("dongjunchung/dpeak")

# ----------------------------------------------------------------
# RNA-seq analysis
# ----------------------------------------------------------------

#cufflinks:
#	cd $(SOURCE_DIR) && \
#	wget -nc http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-$(CUFFLINKS_VER).Linux_x86_64.tar.gz
#	tar xvzf cufflinks-$(CUFFLINKS_VER).Linux_x86_64.tar.gz
#	cp `find $(SOURCE_DIR)/cufflinks-$(CUFFLINKS_VER).Linux_x86_64/ -maxdepth 1 -executable -type f` $(BIN_DIR) && \


ngs_tools: \
	deeptools\
	samtools\
	bedtools\
	sratoolkit\
	fastqc\
	sickle\
	bowtie\
	bowtie2\
	bwa\
	subread\
	tophat2\
	macs1\
	macs2\
	homer

# BEDOPS, DEEPTOOLS
# ================================================================
# Visualization
# ================================================================

# switch to java 8? I recall having an issue with v9, maybe in IGV?
java:
	echo debconf shared/accepted-oracle-license-v1-1 select true | sudo debconf-set-selections
	echo debconf shared/accepted-oracle-license-v1-1 seen true | sudo debconf-set-selections
	sudo add-apt-repository -y ppa:webupd8team/java
	sudo apt-get update
	sudo apt-get -y install oracle-java$(JAVA_VER)-installer
	sudo apt-get install oracle-java$(JAVA_VER)-set-default


igv:
	cd $(SOURCE_DIR) && \
	wget --no-clobber http://data.broadinstitute.org/igv/projects/downloads/IGV_$(IGV_VER).zip && \
	unzip IGV_$(IGV_VER).zip && \
	ln -s -f $(SOURCE_DIR)/IGV_$(IGV_VER)/igv.sh $(BIN_DIR)/igv

igv_tools:
	cd $(SOURCE_DIR) && \
	wget --no-clobber http://data.broadinstitute.org/igv/projects/downloads/igvtools_$(IGVTOOLS_VER).zip  && \
	unzip igvtools_$(IGVTOOLS_VER).zip && \
	ln -s -f $(SOURCE_DIR)/IGVTools/igvtools $(BIN_DIR)/igvtools \
	ln -s -f $(SOURCE_DIR)/IGVTools/igvtools.jar $(BIN_DIR)/igvtools.jar

visualization: java igv igv_tools


## to be used only in non-graphical environment (eg IFB cloud VMs)
desktop_and_x2go:
	sudo apt-get install -y x2goserver
	sudo apt-get install -y --no-install-recommends ubuntu-mate-core ubuntu-mate-desktop
	sudo apt-get install -y mate-desktop-environment-extra
	sudo apt-get install -y mate-notification-daemon caja-gksu caja-open-terminal
	sudo apt-get install -y ambiance-colors radiance-colors


# ================================================================
# All targets
# ================================================================

sequanix_demo: \
	create_bin \
	add_packages \
	R_installation \
	python \
	fastqc \
	bedtools \
	samtools \
	bowtie2 \
	macs2 \
	homer


all: \
	init\
	ubuntu_packages\
	R\
	python\
	ngs_tools\
	visualization

