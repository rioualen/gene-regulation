Dependencies
================================================================

These manuals aim at helping you install the necessary programs and
dependencies in order to have the snakemake workflows work. It was
designed for Unix-running computers (Ubuntu, Debian).

Manual installation
----------------------------------------------------------------

The following manual is meant to help you install the programs that you might need in order to run workflows. 

General requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generic tools
****************************************************************


rsync
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`rsync <https://rsync.samba.org/>`__ is an open source utility that
provides fast incremental file transfer.

::

    sudo apt-get install rsync

git
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Install git on your machine.

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

Several tools require this dependency (e.g. sickle, bamtools...).

::

    sudo apt-get install libz-dev

qsub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



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

Snakemake workflows basic requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python
****************************************************************

Snakemake requires to have Python 3.3+ installed. 
You can check this by issuing the following commands in a terminal:

::

    python --version # usually the default python version is 2.7+
    python3 --version

If you don't have python 3 you should install it.

::

    sudo apt-get install python3

Install pip and pip3.

::

    sudo apt-get install python-pip
    sudo apt-get install python3-pip

Not installed natively?

::

    apt-get install python-dev
    apt-get install python3.4-dev

Pandas library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This library is used in order to read tab-delimited files used in the workflows 
(see files ``samples.tab`` and ``design.tab``).

::

    pip3 install pandas

Package rpy2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    pip3 install "rpy2<2.3.10"



R
****************************************************************


*todo*

Snakemake
****************************************************************

-  `Documentation <http://snakemake.readthedocs.io>`__
-  `FAQ <https://bitbucket.org/snakemake/snakemake/wiki/FAQ>`__
-  `Forum <https://groups.google.com/forum/#!forum/snakemake>`__
-  See also Snakemake section for tutorials. 

Now you have installed Python 3 and pip3 (see previous section), you can
install snakemake safely.

::

    pip3 install snakemake

You can check that snakemake works properly with this basic script:

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

-  Save it to ``~/workspace/hello.py``.
-  Issue the command ``cd ~/workspace ; snakemake -s hello.py``.
-  2 files should be created: ``hello.txt`` and ``bye.txt``.

As of December 2015, you need snakemake version 3.4+.

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

NGS analysis software & tools
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



Alignment/mapping
****************************************************************

The point of mapping is to replace the reads obtained from the sequencing step onto a reference genome. 
When the read is long enough, it can be mapped on the genome with a pretty good confidence, by tolerating a certain amount of so-called mismatches. 
However, genomes can contain repeated regions that are harder to deal with. 

We call "sequencing depth" the average number of reads mapped at each position of the genome. 
The bigger the sequencing depth, the better the quality of the alignment, and the better the downstream analyses. 

BWA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

BWA is designed for short reads alignment. 


`BWA <http://bio-bwa.sourceforge.net/>`__ is a software package for
mapping low-divergent sequences against a large reference genome, such
as the human genome.

-  `Manual <http://bio-bwa.sourceforge.net/bwa.shtml>`__

-  `Publication <http://www.ncbi.nlm.nih.gov/pubmed/19451168>`__ 

Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

::

    sudo apt-get install bwa

Bowtie
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie performs ungapped alignment, and is therefore not suitable for certain types of data, like RNA-seq data. 


::

    cd $HOME/app_sources
    wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip
    unzip bowtie-1.1.1-linux-x86_64.zip
    cp `find bowtie-1.1.1/ -maxdepth 1 -executable -type f` $HOME/bin


Bowtie2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`General
documentation <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__

`Instructions <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2>`__

`Downloads <https://sourceforge.net/projects/bowtie-bio/files/bowtie2/>`__

::

    cd $HOME/app_sources
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip
    unzip bowtie2-2.2.6-linux-x86_64.zip
    p `find bowtie2-2.2.6/ -maxdepth 1 -executable -type f` $HOME/bin


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


