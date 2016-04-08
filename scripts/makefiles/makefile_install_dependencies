export LC_ALL=C
export LANG=C

.PHONY:  \
	usr_bin \
	ngs_bashrc \
	add_pub_key \
	add_repos \
	essential_packages \
	R_installation \
	R_lib \
	python \
	java9 \
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
	homer
	swembl
	igv \
	igv_tools \
	desktop_and_x2go \
	Rstudio 


################################################################
## List targets
MAKEFILE=scripts/makefiles/makefile_install_dependencies
usage:
	echo "usage: make [-OPT='options'] target"
	echo "implemented targets"
	perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

##########################
# Installation variables
##########################

USR_BIN=$(HOME)/bin/
NGS_BASHRC=$(USR_BIN)ngs_bashrc
PATH_NEW=$(PATH)

UBUNTU_VER_NAME="trusty"
CRAN_MIRROR="http://cran.univ-lyon1.fr"
CRAN_PACK_LIST='XML', 'bPeaks', 'caTools'
BIOC_PACK_LIST='affy', 'biomaRt', 'Rsamtools'
PUB_KEY=51716619E084DAB9 F7B8CEA6056E8E56   

###########################
# Tool versions
###########################



BEDTOOLS_VER=2.24.0
BOWTIE1_VER=1.1.1
BOWTIE2_VER=2.2.6
#BWA_VER=0.7.10
FASTQC_VER=0.11.5
IGV_VER=2.3.59
IGVTOOLS_VER=2.3.57
JAVA_VER=oracle-java9
MACS1_VER=1.4.2
NUMPY_VER=1.9.2
RPY2_VER=2.5.6
RSTUDIO_VER=0.99.473
SAMTOOLS_VER=1.3
SRATOOLKIT_VER=2.5.2
SCIPY_VER=0.16.0
SPP_VER=1.11
SWEMBL_VER=3.3.1



##########################
# ~/.bashrc config
##########################


usr_bin:
	mkdir -p $(USR_BIN)
	PATH_NEW=$(PATH_NEW):$(USR_BIN)

ngs_bashrc:
	touch $(NGS_BASHRC)
	echo "" >> $(NGS_BASHRC)
	echo 'alias ls="ls --color"' >> $(NGS_BASHRC)
	echo 'alias rm="rm -i"' >> $(NGS_BASHRC)
	echo "" >> $(NGS_BASHRC)
	echo "export LC_ALL=C" >> $(NGS_BASHRC)
	echo "export LANG=C" >> $(NGS_BASHRC)
	echo 'export force_color_prompt=yes' >> $(NGS_BASHRC)
	echo 'export LC_COLLATE=C' >> $(NGS_BASHRC)
	echo "" >> $(NGS_BASHRC)
#	echo 'export PATH='$(PATH):$(USR_BIN) >> $(NGS_BASHRC)
#	echo 'source $(NGS_BASHRC)' >> /etc/profile



add_pub_key:
	for i in '$(PUB_KEY)'; do echo "PUB_KEY: $$i"; sudo apt-key adv --recv-keys --keyserver keyserver.ubuntu.com $$i; done

add_repos: 
	sudo apt-get install --yes python-software-properties
	sudo add-apt-repository ppa:x2go/stable --yes
	sudo apt-add-repository ppa:ubuntu-mate-dev/ppa --yes
	sudo apt-add-repository ppa:ubuntu-mate-dev/trusty-mate --yes
	sudo add-apt-repository ppa:ravefinity-project/ppa --yes
	sudo apt-get update --yes

essential_packages:
	sudo apt-get update
	sudo apt-get -y install ssh rsync git graphviz gedit-plugins
	sudo apt-get -y install gedit-plugins						
	sudo apt-get -y install zlibc zlib1g-dev				# Required by sickle, bamtools...
	sudo apt-get -y install build-essential						 # Includes gcc compiler
	sudo apt-get -y install libncurses5-dev libncursesw5-dev		# Required at least by samtools
	sudo apt-get -y install libboost-dev						# might me required for spp, to be checked
	sudo apt-get -y install gdebi								# required by rstudio install

desktop_and_x2go:
	sudo apt-get install -y x2goserver
	sudo apt-get install -y --no-install-recommends ubuntu-mate-core ubuntu-mate-desktop
	sudo apt-get install -y mate-desktop-environment-extra
	sudo apt-get install -y mate-notification-daemon caja-gksu caja-open-terminal
	sudo apt-get install -y ambiance-colors radiance-colors;

R_installation:
#	sudo echo "deb $(CRAN_MIRROR)/bin/linux/ubuntu $(UBUNTU_VER_NAME)/" >> /etc/apt/sources.list
	sudo apt-get update
	sudo apt-get -y install r-base r-base-dev libcurl4-openssl-dev libxml2-dev
	echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" >> ~/.Rprofile


R_lib: 
	sudo Rscript -e "pack.list <- c($(CRAN_PACK_LIST)); \
	pack <- pack.list[!(pack.list %in% installed.packages()[,'Package'])]; \
	if(length(pack)) install.packages(pack); \
	source('http://bioconductor.org/biocLite.R'); \
	pack.list <- c($(BIOC_PACK_LIST)); \
	pack <- pack.list[!(pack.list %in% installed.packages()[,'Package'])]; \
	if(length(pack)) biocLite(pack)"

Rstudio: 
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


java9:
	echo debconf shared/accepted-oracle-license-v1-1 select true | sudo debconf-set-selections
	echo debconf shared/accepted-oracle-license-v1-1 seen true | sudo debconf-set-selections
	sudo add-apt-repository -y ppa:webupd8team/java
	sudo apt-get update
	sudo apt-get -y install oracle-java9-installer
	sudo apt-get install oracle-java9-set-default


# ================================================================
#		  NGS
# ================================================================

# ----------------------------------------------------------------
# File management tools
# ----------------------------------------------------------------

samtools:
	cd $(USR_BIN);\
	wget --no-clobber http://sourceforge.net/projects/samtools/files/samtools/$(SAMTOOLS_VER)/samtools-$(SAMTOOLS_VER).tar.bz2;\
	bunzip2 -f samtools-$(SAMTOOLS_VER).tar.bz2;\
	tar xvf samtools-$(SAMTOOLS_VER).tar;\
	rm samtools-$(SAMTOOLS_VER).tar;\
	cd samtools-$(SAMTOOLS_VER); \
	make ;\
	sudo make install;\
	rm -rf samtools*

bedtools:
	cd $(USR_BIN);\
	wget --no-clobber https://github.com/arq5x/bedtools2/releases/download/v$(BEDTOOLS_VER)/bedtools-$(BEDTOOLS_VER).tar.gz;\
	tar xvfz bedtools-$(BEDTOOLS_VER).tar.gz;\
	rm bedtools-$(BEDTOOLS_VER).tar.gz;\
	cd bedtools2; \
	make ;\
	sudo make install;\
	rm -rf bedtools*

sratoolkit:
	wget -nc http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$(SRATOOLKIT_VER)/sratoolkit.$(SRATOOLKIT_VER)-ubuntu64.tar.gz
	mv sratoolkit.$(SRATOOLKIT_VER)-ubuntu64.tar.gz $(USR_BIN)
	cd $(USR_BIN); \
	tar xvzf sratoolkit.$(SRATOOLKIT_VER)-ubuntu64.tar.gz; \
	rm sratoolkit.$(SRATOOLKIT_VER)-ubuntu64.tar.gz; \
	ln -s -f sratoolkit.$(SRATOOLKIT_VER)-ubuntu64/bin/fastq-dump fastq-dump; \
	ln -s -f sratoolkit.$(SRATOOLKIT_VER)-ubuntu64/bin/prefetch prefetch


# ----------------------------------------------------------------
# Quality assessment & trimming
# ----------------------------------------------------------------

fastqc:
#	sudo apt-get -y install fastqc
	wget --no-clobber http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v$(FASTQC_VER).zip;\
	mv fastqc_v$(FASTQC_VER).zip $(USR_BIN); \
	cd $(USR_BIN); \
	unzip -o fastqc_v$(FASTQC_VER).zip; \
	rm fastqc_v$(FASTQC_VER).zip; \
	sudo chmod +x FastQC/fastqc; \
	ln -s -f FastQC/fastqc fastqc; \

sickle: 
	git clone https://github.com/najoshi/sickle.git; \
	cd sickle;\
	make; \
	mv sickle $(USR_BIN) ;\
	cd ..;\
	rm -Rf sickle

# ----------------------------------------------------------------
# Mapping tools
# ----------------------------------------------------------------

bowtie: 
	cd $(USR_BIN); \
	wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie/$(BOWTIE1_VER)/bowtie-$(BOWTIE1_VER)-linux-x86_64.zip;\
	unzip bowtie-$(BOWTIE1_VER)-linux-x86_64.zip;\
	rm bowtie-$(BOWTIE1_VER)-linux-x86_64.zip; \
#	echo 'export PATH='$(PATH):$(USR_BIN)bowtie-$(BOWTIE1_VER) >> $(NGS_BASHRC)
	PATH_NEW=$(PATH_NEW):$(USR_BIN)bowtie-$(BOWTIE1_VER);\

bowtie2:
	cd $(USR_BIN); \
	wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/$(BOWTIE2_VER)/bowtie2-$(BOWTIE2_VER)-linux-x86_64.zip;\
	unzip bowtie2-$(BOWTIE2_VER)-linux-x86_64.zip;\
	rm bowtie2-$(BOWTIE2_VER)-linux-x86_64.zip;\
#	echo 'export PATH='$(PATH):$(USR_BIN)bowtie2-$(BOWTIE2_VER) >> $(NGS_BASHRC)
	PATH_NEW=$(PATH_NEW):$(USR_BIN)bowtie2-$(BOWTIE2_VER);\

bwa:
	sudo apt-get -y install bwa
#	wget -nc https://sourceforge.net/projects/bio-bwa/files/bwa-$(BWA_VER).tar.bz2; \




# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------

macs1:
	cd $(USR_BIN); \
	wget --no-clobber https://github.com/downloads/taoliu/MACS/MACS-$(MACS1_VER)-1.tar.gz; \
	tar xvfz MACS-$(MACS1_VER)-1.tar.gz ; \
	rm MACS-$(MACS1_VER)-1.tar.gz ; \
	cd MACS-$(MACS1_VER); \
	sudo python setup.py install

macs2:
	sudo pip install MACS2;

spp:
	cd $(USR_BIN);\
	wget http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_$(SPP_VER).tar.gz;\
	sudo R CMD INSTALL spp_$(SPP_VER).tar.gz; \

swembl:
	cd $(USR_BIN);\
	wget "http://www.ebi.ac.uk/~swilder/SWEMBL/SWEMBL.$(SWEMBL_VER).tar.bz2"; \
	bunzip2 -f SWEMBL.$(SWEMBL_VER).tar.bz2;\
	tar xvf SWEMBL.$(SWEMBL_VER).tar;\
	rm SWEMBL.$(SWEMBL_VER).tar;\
	chown -R ubuntu-user SWEMBL.$(SWEMBL_VER);\
	cd SWEMBL.$(SWEMBL_VER);\
	make

homer:
	wget "http://homer.salk.edu/homer/configureHomer.pl"; \
	mkdir $(USR_BIN)/HOMER; \
	mv configureHomer.pl $(USR_BIN)/HOMER; \
	cd $(USR_BIN)/HOMER; \
	perl configureHomer.pl -install homer; \
#	echo 'export PATH='$(PATH):$(USR_BIN)HOMER/bin >> $(NGS_BASHRC)
	PATH_NEW=$(PATH_NEW):$(USR_BIN)HOMER/bin



# ================================================================
#		  Visualization
# ================================================================

igv:
	cd $(USR_BIN); \
	wget --no-clobber http://data.broadinstitute.org/igv/projects/downloads/IGV_$(IGV_VER).zip; \
	unzip IGV_$(IGV_VER).zip;\
	rm  IGV_$(IGV_VER).zip;\
	ln -s -f $(USR_BIN)/IGV_$(IGV_VER)/igv.sh $(USR_BIN)/igv

igv_tools:
	wget --no-clobber http://data.broadinstitute.org/igv/projects/downloads/igvtools_$(IGVTOOLS_VER).zip ;\
	mv igvtools_$(IGVTOOLS_VER).zip $(USR_BIN);\
	cd $(USR_BIN); \
	unzip igvtools_$(IGVTOOLS_VER).zip;\
	rm igvtools_$(IGVTOOLS_VER).zip;\
	ln -s -f $(USR_BIN)/IGVTools/igvtools $(USR_BIN)/igvtools

edit_path:
	echo 'export PATH='$(PATH_NEW) >> $(NGS_BASHRC)

all: \
	usr_bin \
	ngs_bashrc \
	add_pub_key \
	add_repos \
	essential_packages \
	R_installation \
	R_lib \
	python \
	java9 \
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
	igv \
	igv_tools \
	desktop_and_x2go \
	edit_path

