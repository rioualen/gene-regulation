Dependencies
================================================================

These manuals aim at helping you install programs and
dependencies used in the Gene-regulation library. 

Some of them are mandatory, and some are optional, depending 
on the Snakemake workflows you need to run. 

They were tested under Ubuntu 14.04. 

Manual installation
----------------------------------------------------------------

This manual is organized in sections, so you can cherry-pick the programs you want to manually install. 
For "all inclusive" solutions, please refer yourself to the following sections. 

General 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generic tools
****************************************************************

nano
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nano is a simple command-line text editor. 

::

    sudo apt-get install nano


rsync
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Rsync <https://rsync.samba.org/>`__ is an open source utility that
provides fast incremental file transfer.

::

    sudo apt-get install rsync

git
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Git rsync <https://en.wikipedia.org/wiki/Git>`__ is a version control system (VCS) for tracking changes in computer files and coordinating work on those files among multiple people. 

::

    sudo apt-get install git

Optional:

-  Create an account on `GitHub <https://github.com>`__.
-  Add your ssh public key to your GitHub account settings (account >
   settings > SSH keys > add SSH key).

::

    less ~/.ssh/id_rsa.pub

zlib
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unix package required by several tools, including Sickle and Bamtools.

::

    sudo apt-get install libz-dev

Java
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Java is required by several tools using GUIs, such as FastQC or IGV. 

It seems java 9 causes issues with IGV, so we chose to use java 8 here. 

::

	echo debconf shared/accepted-oracle-license-v1-1 select true | sudo debconf-set-selections
	echo debconf shared/accepted-oracle-license-v1-1 seen true | sudo debconf-set-selections
	sudo add-apt-repository -y ppa:webupd8team/java
	sudo apt-get update
	sudo apt-get -y install oracle-java8-installer

Check installation:

::

     java -version

Create bin/ and app\_sources/ (optional)
****************************************************************

While some programs will be installed completely automatically, others 
will not. Here we create a directory that will be used for manual
installations.

::

    mkdir $HOME/bin
    mkdir $HOME/app_sources

You might then have to edit your ``$PATH`` manually (see next section).

Edit ``$PATH``
****************************************************************

In order to use manually installed programs and make them executable,
you may have to update your ``$PATH`` environment variable. You can do
so by editing the ``~/.profile`` file.

::

    nano ~/.profile

Fetch this paragraph and add the path to manually installed executables:

::

    # set PATH so it includes user's private bin if it exists
    if [ -d "$HOME/bin" ] ; then
        PATH="$HOME/bin:$PATH"
    fi

Execute the file to validate the change.

::

    source ~/.profile

Snakemake basic requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python
****************************************************************

Snakemake requires to have Python 3.3+ installed. 
You can check this by issuing the following commands in a terminal:

::

    python --version
    python3 --version

If you don't have python 3 you should install it.

::

    sudo apt-get install python3

Install python package managers and devel libraries.

::

    apt-get install python-dev
    apt-get install python3.4-dev
    sudo apt-get install python-pip
    sudo apt-get install python3-pip


Pandas library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Python Data Analysis Library <http://pandas.pydata.org/>`__ is an open source, BSD-licensed library providing high-performance, easy-to-use data structures and data analysis tools for the Python programming language.

This library is used in order to read tab-delimited files used in the workflows 
(see files ``samples.tab`` and ``design.tab``).

::

    pip3 install pandas

Package rpy2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The package `rpy2 <https://rpy2.readthedocs.io>`__ alloàws to access R from within Python code. 

::

    pip3 install "rpy2<2.3.10"



R
****************************************************************

You can fetch a CRAN mirror `here <https://cran.r-project.org/mirrors.html>`__. 

::

	sudo sh -c "echo 'deb <your mirror> trusty/' >> /etc/apt/sources.list"                          ## Repository for Ubuntu 14.04 Trusty Tahr
	#sudo sh -c "echo 'deb http://ftp.igh.cnrs.fr/pub/CRAN/ trusty/' >> /etc/apt/sources.list"      ## Mirror in Montpellier, France
	sudo apt-get -y update
	sudo apt-get -y install r-base r-base-dev libcurl4-openssl-dev libxml2-dev
	echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" >> ~/.Rprofile

Check installation:

::

    R --version

Snakemake
****************************************************************

"Snakemake is a workflow engine that provides a readable Python-based workflow definition language and a powerful execution environment 
that scales from single-core workstations to compute clusters without modifying the workflow. 
It is the first system to support the use of automatically inferred multiple named wildcards (or variables) in input and output filenames."

(Köster and Rahman, 2012)

-  `Documentation <http://snakemake.readthedocs.io>`__
-  `FAQ <https://bitbucket.org/snakemake/snakemake/wiki/FAQ>`__
-  `Forum <https://groups.google.com/forum/#!forum/snakemake>`__
-  See also Snakemake section for tutorials. 

NB: Python 3 and pip3 are required ('see `this section <http://gene-regulation.readthedocs.io/en/latest/dependencies.html#python>`__). 

::

    pip3 install snakemake

You can check that snakemake works properly with this basic script. 

::

    """Snakefile to test basic functions of snakemake.
    """
    rule all:
        input: expand("bye.txt")

    rule hello:
        """Write HELLO in a text file named hello.txt.
        """
        output: "hello.txt"
        message: "Generating {output} file."
        shell: "echo HELLO > {output}"

    rule bye:
        """Write BYE in a text file named bye.txt.
        """
        input: "hello.txt"
        output: "bye.txt"
        message: "Generating {output} file."
        shell: "echo BYE > {output}"

::

    touch $HOME/hello.py
    nano $HOME/hello.py             ## copy/paste script above and save

Execute the workflow; two files should be created: ``hello.txt`` and ``bye.txt``.

::

    cd ; snakemake -s hello.py

In case you need to upgrade snakemake:

::

    pip3 install snakemake --upgrade

If you want to use Snakemake reports function (optional):

::

    pip3 install docutils

Graphviz
****************************************************************

Snakemake can generate useful graphviz outputs.

::

    sudo apt-get install graphviz

NGS analysis tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quality assessment
****************************************************************

FastQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__
aims to provide a simple way to do some quality control checks on raw
sequence data coming from high throughput sequencing pipelines. It
provides a modular set of analyses which you can use to give a quick
impression of whether your data has any problems of which you should be
aware before doing any further analysis.

The main functions of FastQC are:

-  Import of data from BAM, SAM or FastQ files (any variant)
-  Providing a quick overview to tell you in which areas there may be
   problems
-  Summary graphs and tables to quickly assess your data
-  Export of results to an HTML based permanent report
-  Offline operation to allow automated generation of reports without
   running the interactive application

Links:

-  `QC Fail Sequencing <https://sequencing.qcfail.com/>`__

-  `FastQC results
   interpretation <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/>`__

FastQC is available from linux repositories:

::

    sudo apt-get install fastqc

However, since it's an older version, it can cause problems of dependencies. 

We recommend installing it manually: 

::

    cd $HOME/app_sources
    wget --no-clobber http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
    unzip -o fastqc_v0.11.5.zip
    chmod +x FastQC/fastqc
    ln -s -f $HOME/app_sources/FastQC/fastqc $HOME/bin/fastqc

NB: FastQC requires to have Java installed (even for commandline use). 

Check installation:

::

    fastqc --version

MultiQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`MultiQC <http://multiqc.info/>`__ searches a given directory for analysis logs and compiles a HTML report. 
It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.

::

    sudo pip install multiqc

NB: a bug can appear depending on versions:

::

Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/matplotlib
Storing debug log for failure in /home/gr/.pip/pip.log

If so, it can be avoided by installing ubuntu dependencies, then reinstalling multiqc:

::

    sudo apt-get install libfreetype6-dev python-matplotlib
    sudo pip install multiqc

Check installation: 

::

     multiqc --version

Trimming
****************************************************************

The quality of the reads generated by high-throughput sequencing
technologies tends to decrease at their ends. Trimming consists in
cutting out theses ends, and thus better the quality of reads before the
mapping.


Sickle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Sickle <https://github.com/najoshi/sickle>`__ is a trimming tool which
better the quality of NGS reads.

Sickle uses sliding windows computing sequencing quality along the
reads. When the quality falls below a chose q-value threshold, the reads
is cut. If the size of the remaining read is too short, it is completely
removed. Sickle takes into account three different types of read
quality: Illumina, Solexa, Sanger.


-  Pre-requisite: install ``zlib`` (*link to section*).
-  Clone the git repository into your bin (*link to section*) and run
   ``make``.

::

    cd $HOME/app_sources
    git clone https://github.com/najoshi/sickle.git 
    cd sickle 
    make 
    cp sickle $HOME/bin

Check installation: 

::

    sickle --version

Cutadapt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Cutadapt <http://cutadapt.readthedocs.io/en/stable/>`__ finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

::

    pip install --user --upgrade cutadapt

Check installation:

::

    cutadapt --version


TrimGalore
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In our workflows we use `TrimGalore <https://github.com/FelixKrueger/TrimGalore>`__, a wrapper around Cutadapt and FastQC. 
It should be installed if you want to run cutadapt. 

::

    cutadapt --version                              # Check that cutadapt is installed
    fastqc -v                                       # Check that FastQC is installed

    cd $HOME/app_sources
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz -o trim_galore.tar.gz
    tar xvzf trim_galore.tar.gz
    mv TrimGalore-0.4.3/trim_galore $HOME/bin

Check installation:

::

    trim_galore --version


Alignment/mapping
****************************************************************

The point of mapping is to replace the reads obtained from the sequencing step onto a reference genome. 
When the read is long enough, it can be mapped on the genome with a pretty good confidence, by tolerating a certain amount of so-called mismatches. 
However, genomes can contain repeated regions that are harder to deal with. 

We call "sequencing depth" the average number of reads mapped at each position of the genome. 
The bigger the sequencing depth, the better the quality of the alignment, and the better the downstream analyses. 

BWA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`BWA <http://bio-bwa.sourceforge.net/>`__ is a software package for
mapping low-divergent sequences against a large reference genome, such
as the human genome. It is designed for short reads alignment. 


-  `Manual <http://bio-bwa.sourceforge.net/bwa.shtml>`__

-  `Publication <http://www.ncbi.nlm.nih.gov/pubmed/19451168>`__ 

Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

::

    sudo apt-get install bwa

Check installation:

::

    bwa

Bowtie
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie performs ungapped alignment, and is therefore not suitable for certain types of data, like RNA-seq data. 


::

    cd $HOME/app_sources
    wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip
    unzip bowtie-1.1.1-linux-x86_64.zip
    cp `find bowtie-1.1.1/ -maxdepth 1 -executable -type f` $HOME/bin

Check installation:

::

     bowtie --help

Bowtie2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"`Bowtie 2 <http://bowtie-bio.sourceforge.net>`__ is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. 
It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters to relatively long (e.g. mammalian) genomes. 
Bowtie 2 indexes the genome with an FM Index (based on the Burrows-Wheeler Transform or BWT) to keep its memory footprint small: 
for the human genome, its memory footprint is typically around 3.2 gigabytes of RAM. 
Bowtie 2 supports gapped, local, and paired-end alignment modes. 
Multiple processors can be used simultaneously to achieve greater alignment speed. 
Bowtie 2 outputs alignments in SAM format, enabling interoperation with a large number of other tools (e.g. SAMtools, GATK) that use SAM. 
Bowtie 2 is distributed under the GPLv3 license, and it runs on the command line under Windows, Mac OS X and Linux."

`General
documentation <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__

`Instructions <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2>`__

`Downloads <https://sourceforge.net/projects/bowtie-bio/files/bowtie2/>`__

Reference:

Langmead B, Trapnell C, Pop M, L Salzberg S. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 200910:R25. DOI: 10.1186/gb-2009-10-3-r25


::

    cd $HOME/app_sources
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip
    unzip bowtie2-2.2.6-linux-x86_64.zip
    cp `find bowtie2-2.2.6/ -maxdepth 1 -executable -type f` $HOME/bin

Subread-align
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
http://subread.sourceforge.net/


he Subread package comprises a suite of software programs for processing next-gen sequencing read data including:

    Subread: a general-purpose read aligner which can align both genomic DNA-seq and RNA-seq reads. It can also be used to discover genomic mutations including short indels and structural variants.
    Subjunc: a read aligner developed for aligning RNA-seq reads and for the detection of exon-exon junctions. Gene fusion events can be detected as well.
    featureCounts: a software program developed for counting reads to genomic features such as genes, exons, promoters and genomic bins.
    exactSNP: a SNP caller that discovers SNPs by testing signals against local background noises

Ref

Liao Y, Smyth GK and Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013
Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014


subread:
	cd $(SOURCE_DIR) && \
	wget -nc https://sourceforge.net/projects/subread/files/subread-$(SUBREAD_VER)/subread-$(SUBREAD_VER)-source.tar.gz && \
	tar zxvf subread-$(SUBREAD_VER)-source.tar.gz && \
	cd subread-$(SUBREAD_VER)-source/src && \
	make -f Makefile.Linux && \
	cd ../bin && \
	cp `find * -executable -type f` $(BIN_DIR)




Tophat
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Peak-calling
****************************************************************

The following tools can be used to perform ChIP-seq peak-calling.


HOMER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Required in order to run the demo workflow "ChIP-seq" on dataset GSE20870 (in the tutorials section). 

`Web page <http://homer.salk.edu/>`__

`Install
instructions <http://homer.salk.edu/homer/introduction/install.html>`__

::

    mkdir $HOME/app_sources/homer
    cd $HOME/app_sources/homer
    wget "http://homer.salk.edu/homer/configureHomer.pl"
    perl configureHomer.pl -install homer
    cp `find $HOME/app_sources/homer/bin -maxdepth 1 -executable -type f` $HOME/bin

The basic Homer installation does not contain any sequence data. To
download sequences for use with HOMER, use the configureHomer.pl script.
To get a list of available packages:

::

    perl $HOME/bin/HOMER/configureHomer.pl -list

To install packages, simply use the -install option and the name(s) of
the package(s).

::

    perl  $HOME/bin/HOMER/configureHomer.pl -install mouse # (to download the mouse promoter set)
    perl  $HOME/bin/HOMER/configureHomer.pl -install mm8   # (to download the mm8 version of the mouse genome)
    perl  $HOME/bin/HOMER/configureHomer.pl -install hg19  # (to download the hg19 version of the human genome)

Supported organisms:

+-----------------+--------------------+
| Organism        | Assembly           |
+=================+====================+
| Human           | hg17, hg18, hg19   |
+-----------------+--------------------+
| Mouse           | mm8, mm9, mm10     |
+-----------------+--------------------+
| Rat             | rn4, rn5           |
+-----------------+--------------------+
| Frog            | xenTro2, xenTro3   |
+-----------------+--------------------+
| Zebrafish       | danRer7            |
+-----------------+--------------------+
| Drosophila      | dm3                |
+-----------------+--------------------+
| C. elegans      | ce6, ce10          |
+-----------------+--------------------+
| S. cerevisiae   | sacCer2, sacCer3   |
+-----------------+--------------------+
| S. pombe        | ASM294v1           |
+-----------------+--------------------+
| Arabidopsis     | tair10             |
+-----------------+--------------------+
| Rice            | msu6               |
+-----------------+--------------------+

HOMER can also work with custom genomes in FASTA format and gene
annotations in GTF format.

MACS 1.4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Required in order to run the demo workflow "ChIP-seq" on dataset GSE20870 (in the tutorials section). 


-  `Documentation <http://liulab.dfci.harvard.edu/MACS/00README.html>`__
-  `Installation manual <http://liulab.dfci.harvard.edu/MACS/INSTALL.html>`__

::

    cd $HOME/app_sources
    wget "https://github.com/downloads/taoliu/MACS/MACS-1.4.3.tar.gz"
    tar -xvzf MACS-1.4.3.tar.gz
    cd MACS-1.4.3
    sudo python setup.py install
    macs14 --version


MACS2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Required in order to run the demo workflow "ChIP-seq" on dataset GSE20870 (in the tutorials section). 

-  `Webpage <https://github.com/taoliu/MACS/>`__

::

    sudo apt-get install python-numpy
    sudo pip install MACS2

bPeaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-caller developped specifically for yeast, can be useful in order to
process small genomes only.

It is currently not used in demo workflows, and is therefore not m adatory to run the tutorials. 

Available as an R library.

`Web page <http://bpeaks.gene-networks.net/>`__

::

    install.packages("bPeaks")
    library(bPeaks)


SPP R package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation of this peak-caller is optional, as it is not currently published and maintained properly. 

It is therefore not used in our demo workflows. 


- In R

::

    source("http://bioconductor.org/biocLite.R")
    biocLite("spp")
    install.packages("caTools")
    install.packages("spp")

- In commandline

::

    apt-get install libboost-all-dev
    cd $HOME/app_sources
    wget -nc http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.11.tar.gz
    sudo R CMD INSTALL spp_1.11.tar.gz

- Using git (I haven't tried this one but it looks more recent) (see `github page <https://github.com/hms-dbmi/spp>`__)

::

    require(devtools)
    devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)


I also wrote a little protocol a while ago. 
Here's the procedure on Ubuntu 14.04, in this very order:

In unix shell:

::

    # unix libraries
    apt-get update
    apt-get -y install r-base
    apt-get -y install libboost-dev zlibc zlib1g-dev

In R shell:

::

    # Rsamtools
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")

In unix shell:

::

    # spp
    wget http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.11.tar.gz
    sudo R CMD INSTALL spp_1.11.tar.gz

A few links:

-  Download page can be found
   `here <http://compbio.med.harvard.edu/Supplements/ChIP-seq/>`__,
   better chose version ``1.11``.
-  SPP requires the Bioconductor library
   `Rsamtools <https://bioconductor.org/packages/release/bioc/html/Rsamtools.html>`__
   to be installed beforehand.
-  Unix packages ``gcc`` and ``libboost`` (or equivalents) must be
   installed.
-  You can find a few more notes
   `here <http://seqanswers.com/forums/archive/index.php/t-22653.html>`__.
-  Good luck!

SWEMBL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation of this peak-caller is optional, as it is not currently published and maintained properly. 

It is therefore not used in our demo workflows. 


-  `SWEMBL beginner's
   manual <http://www.ebi.ac.uk/~swilder/SWEMBL/beginners.html>`__

::

    cd $HOME/app_sources
    wget "http://www.ebi.ac.uk/~swilder/SWEMBL/SWEMBL.3.3.1.tar.bz2"
    bunzip2 -f SWEMBL.3.3.1.tar.bz2
    tar xvf SWEMBL.3.3.1.tar
    rm SWEMBL.3.3.1.tar
    chown -R ubuntu-user SWEMBL.3.3.1
    cd SWEMBL.3.3.1
    make

It seems there could be issues with C flags. To be investigated. 

Motif discovery, motif analysis
****************************************************************

These software can be useful for the analysis of ChIP-seq peaks. 

Regulatory Sequence Analysis Tools (RSAT)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*see dedicated section*

`Link <http://rsat.eu/>`__

*to translate*

Suite logicielle spécialisée pour l'analyse de motifs cis-régulateurs,
développée par les équipes de Morgane Thomas-Chollier (ENS, Paris) et
Jacques van Helden (TAGC, Marseille). Inclut des outils spécifiques pour
l'analyse de données de ChIP-seq.



MEME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Link <http://meme.ebi.edu.au/meme/doc/meme-chip.html>`__

*to translate*

Suite logicielle spécialisée pour l'analyse de motifs cis-régulateurs,
développée par l'équipe de Tim Bailey. Inclut des outils spécifiques
pour l'analyse de données de ChIP-seq.


Miscellaneous
****************************************************************

SRA toolkit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This toolkit includes a number of programs, allowing the conversion of
``*.sra`` files. ``fastq-dump`` translates ``*.sra`` files to
``*.fastq`` files.

-  `SRA format <http://www.ncbi.nlm.nih.gov/Traces/sra/>`__
-  `fastq-dump
   manual <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump>`__
-  `Installation
   manual <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std>`__

You can download last version
`here <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__,
or issue the following commands:

::

    cd $HOME/app_sources
    wget -nc http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz
    tar xzf sratoolkit.2.5.2-ubuntu64.tar.gz
    cp `find sratoolkit.2.5.2-ubuntu64/bin -maxdepth 1 -executable -type l` $HOME/bin

You can also install SRA toolkit simply by issuing this
command, but likely it won't be the most recent release:

::

    sudo apt-get install sra-toolkit

::

    fastq-dump --version
      fastq-dump : 2.1.7

Samtools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SAM (Sequence Alignment/Map) format is a generic format for storing
large nucleotide sequence alignments.

`SAMtools <http://samtools.sourceforge.net/>`__ provides several tools
to process such files.

::

    cd $HOME/app_sources
    wget --no-clobber http://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
    bunzip2 -f samtools-1.3.tar.bz2
    tar xvf samtools-1.3.tar
    cd samtools-1.3
    make 
    sudo make install

Bedtools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `bedtools <http://bedtools.readthedocs.org/en/latest/>`__ utilities
are a swiss-army knife of tools for a wide-range of genomics analysis
tasks. For example, bedtools allows one to intersect, merge, count,
complement, and shuffle genomic intervals from multiple files in
widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF.

::

    sudo apt-get install bedtools

or get the latest version:

::

    cd $HOME/app_sources
    wget --no-clobber https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz
    tar xvfz bedtools-2.24.0.tar.gz
    cd bedtools2
    make
    sudo make install



Bedops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cd $HOME/app_sources
    wget -nc https://github.com/bedops/bedops/releases/download/v2.4.19/bedops_linux_x86_64-v2.4.19.tar.bz2
    tar jxvf bedops_linux_x86_64-v2.4.19.tar.bz2
    mkdir bedops
    mv bin bedops
    cp bedops/bin/* $HOME/bin

Deeptools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cd $HOME/app_sources
    git clone https://github.com/fidelram/deepTools
    cd deepTools
    python setup.py install

Picard tools 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*todo*

Other
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  `MICSA <http://bioinfo-out.curie.fr/software.html>`__: peak-calling &
   motifs discovery
   (`publication <http://bioinformatics.oxfordjournals.org/content/26/20/2622.long>`__).
-  `ChIPMunk <http://line.imb.ac.ru/ChIPMunk>`__: deep and wide digging
   for binding motifs in ChIP-Seq data
   (`publication <http://bioinformatics.oxfordjournals.org/content/26/20/2622.long>`__).
-  `HMCan <http://www.cbrc.kaust.edu.sa/hmcan/>`__: a method for
   detecting chromatin modifications in cancer samples using ChIP-seq
   data
   (`publication <http://bioinformatics.oxfordjournals.org/content/29/23/2979.long>`__).
-  seqMINER
-  `Crunch project <http://crunch.unibas.ch/fcgi/crunch.fcgi>`__
-  CSDeconv
-  ...



Makefile
----------------------------------------------------------------

The Gene-regulation library comprises a makefile that can install most of the 
dependencies described in the previous section. It is recommended when you're setting up a virtual environments, 
as described in `these tutorials <http://gene-regulation.readthedocs.io/en/latest/environments.html>`_. 

If you want to run the workflows on your personal computer or on a server, you should follow the `manual installation 
<http://gene-regulation.readthedocs.io/en/latest/dependencies.html#manual-installation>`_, or contact a sysadmin. 

The makefile currently allows running the following workflows:

- import_from_sra.wf
- quality_control.wf
- ChIP-seq.wf

It is not yet handling al the RNA-seq dependencies.

::

    # it is assumed that you have defined a global variable with the path to the Gene-regulation library
    cd ${GENE_REG_PATH}
    make -f scripts/makefiles/install_tools_and_libs.mk all
    source ~/.bashrc

..
Conda
----------------------------------------------------------------
A number of dependencies of Gene-regulation can be installed through a Conda environment. 
This list is not exhaustive. 
    conda install -c bioconda sickle=0.5 
    conda install -c bioconda bowtie=1.2.0 
    conda install -c bioconda bowtie2=2.3.0 
    conda install -c bioconda subread=1.5.0.post3 
    conda install -c bioconda tophat=2.1.1 
    conda install -c bioconda bwa=0.7.15 
    conda install -c bioconda fastqc=0.11.5 
    conda install -c bioconda macs2=2.1.1.20160309 
    conda install -c bioconda homer=4.8.3 
    conda install -c bioconda bedtools=2.26.0 
    conda install -c bioconda samtools=1.3.1 
    conda install -c bioconda bamtools=2.4.0 


