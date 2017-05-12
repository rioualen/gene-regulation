
This manual aims at helping you install the necessary programs and
dependencies in order to have the snakemake workflows work. It was
designed for Unix-running computers (Ubuntu, Debian).

General requirements
----------------------------------------------------------------

Generic tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ssh
****************************************************************

Secure Shell (SSH) is a network protocol for operating network services securely over an unsecured network.
 The best known example application is for remote login to computer systems by users.

It can be useful, for example, in order to connect to a virtual machine from a host terminal, 
work on a remote server, or connect with a git repository. 


::

    sudo apt-get install ssh

In order to use this protocol, you should dispose of an ssh key. It is usually located in your home directory, in a hidden directory: `$HOME/.ssh`. 
If you don't have a key, and only then, you should create one. Follow for example `this tutorial <https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/>`_. 

rsync
****************************************************************

`rsync <https://rsync.samba.org/>`_ is an open source utility that
provides fast incremental file transfer.

::

    sudo apt-get install rsync

git
****************************************************************

Git is a version control system (VCS) for tracking changes in computer files and coordinating work on those files among multiple people. 
It is primarily used for software development, but it can be used to keep track of changes in any files. 
As a distributed revision control system it is aimed at speed, data integrity, and support for distributed, non-linear workflows. 
(from `Wikipedia page <https://en.wikipedia.org/wiki/Git>`_

::

    sudo apt-get install git

Gene-regulation is distibuted on the `GitHub <https://en.wikipedia.org/wiki/GitHub>`_ hosting service.  

-  Create an account on `GitHub <https://github.com>`_.
-  Install git on your machine.

-  Add your ssh public key to your GitHub account settings (account >
   settings > SSH keys > add SSH key).

::

    less ~/.ssh/id_rsa.pub

zlib
****************************************************************

Several tools require this dependency (e.g. sickle, bamtools...).

::

    sudo apt-get install libz-dev

qsub
****************************************************************

**TODO**

Create bin/ and app\_sources/ (opt) TODO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While some programs will be installed completely automatically... some
will not. Here we create a directory that will be used for manual
installations.

::

    mkdir $HOME/bin
    mkdir $HOME/app_sources

You might then have to edit your ``$PATH`` manually (see next section).

Edit ``$PATH``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to use manually installed programs and make them executable,
you may have to update your ``$PATH`` environment variable. You can do
so by editing the ``~/.profile`` file.

::

    nano ~/.profile

Fetch this paragraph and add your new tool:

::

    # set PATH so it includes user's private bin if it exists
    if [ -d "$HOME/bin" ] ; then
        PATH="$HOME/bin/sickle:$PATH"
    fi

Execute the file to validate the change.

::

    source ~/.profile

Snakemake workflows basic requirements
----------------------------------------------------------------

.. raw:: html

   <!--
   /!\\ check R versions, python versions, etc...

   /!\\ It seems there could an issue with the installation of python module while using sudo or doing is as root...
   Googling "sudo pip" suggests it is a bad idea for safety reasons, and should not be needed...
   -->

**TODO:**

-  R section (including lib path issue)
-  Is running pip/pip3 as sudoer necessary or desirable?

Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snakemake requires to have Python 3.3+ installed. If you have installed
Ubuntu14.04 on your machine by following one of the previous tutorials,
you should have both Python 2.7 and Python 3.4 installed for starters.
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
****************************************************************

This library is used in order to read our data, specifically files
``samples.tab`` and ``design.tab``.

::

    pip3 install pandas

Package rpy2
****************************************************************

::

    pip3 install "rpy2<2.3.10"

**NB** There might be other dependencies; this should be checked by
running workflows without the RSAT install, which itself contains many
libraries.

R
-

**TODO**

.. raw:: html

   <!-- unnecessary -> use rsat ?
   ### Biostrings (peak length)

   ```
   install.packages("Biostrings", lib="/path/to/my/lib")
   ```
   -->

Snakemake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `Documentation <https://bitbucket.org/snakemake/snakemake/wiki/Documentation>`_
-  `FAQ <https://bitbucket.org/snakemake/snakemake/wiki/FAQ>`_
-  `Forum <https://groups.google.com/forum/#!forum/snakemake>`_
-  More: see
   `wiki/informatics <https://github.com/rioualen/gene-regulation/blob/master/doc/wiki-fg/informatics.md>`_
   section.

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snakemake can generate useful graphviz outputs.

::

    sudo apt-get install graphviz

NGS analysis software & tools
----------------------------------------------------------------

File management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SRA toolkit
****************************************************************

This toolkit includes a number of programs, allowing the conversion of
``*.sra`` files. ``fastq-dump`` translates ``*.sra`` files to
``*.fastq`` files.

-  `SRA format <http://www.ncbi.nlm.nih.gov/Traces/sra/>`_
-  `fastq-dump
   manual <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump>`_
-  `Installation
   manual <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std>`_

You can download last version
`here <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_,
or issue the following commands:

::

    cd ~/bin
    wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz"
    tar -xvzf sratoolkit.2.5.2-ubuntu64.tar.gz
    rm sratoolkit.2.5.2-ubuntu64.tar.gz

Add to path (cf section 1.3):

::

    PATH="$HOME/bin/sratoolkit.2.5.2-ubuntu64/bin:$PATH"

Check version:

::

    fastq-dump --version
      fastq-dump : 2.5.2

You should be able to install SRA toolkit simply by issuing this
command, but likely it won't be the most recent release:

::

    sudo apt-get install sra-toolkit

::

    fastq-dump --version
      fastq-dump : 2.1.7

Samtools
****************************************************************

SAM (Sequence Alignment/Map) format is a generic format for storing
large nucleotide sequence alignments.

`SAMtools <http://samtools.sourceforge.net/>`_ provides several tools
to process such files.

TODO: install samtools from website, not from apt-get repositories.

.. raw:: html

   <!--
   ```
   sudo apt-get install samtools
   ```
   V: 0.1.19
   Latest: 1.2
   -->

Bedtools
****************************************************************

The `bedtools <http://bedtools.readthedocs.org/en/latest/>`_ utilities
are a swiss-army knife of tools for a wide-range of genomics analysis
tasks. For example, bedtools allows one to intersect, merge, count,
complement, and shuffle genomic intervals from multiple files in
widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF.

::

    sudo apt-get install bedtools

V: v2.17.0 Latest: 2.24.0

Quality assessment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FastQC
****************************************************************

`FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
aims to provide a simple way to do some quality control checks on raw
sequence data coming from high throughput sequencing pipelines. It
provides a modular set of analyses which you can use to give a quick
impression of whether your data has any problems of which you should be
aware before doing any further analysis.

::

    sudo apt-get install fastqc

Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sickle
****************************************************************

`Sickle <https://github.com/najoshi/sickle>`_ is a trimming tool which
better the quality of NGS reads.

-  Pre-requisite: install ``zlib`` (see section 1.1.4).
-  Clone the git repository into your bin (see section 1.2) and run
   ``make``.

::

    cd $HOME/bin
    git clone https://github.com/najoshi/sickle.git
    cd sickle
    make

-  Add sickle to your ``$PATH`` (see section 1.3).

::

    PATH="$HOME/bin/sickle:$PATH"

Alignment/mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BWA
****************************************************************

`BWA <http://bio-bwa.sourceforge.net/>`_ is a software package for
mapping low-divergent sequences against a large reference genome, such
as the human genome.

-  `Manual <http://bio-bwa.sourceforge.net/bwa.shtml>`_

::

    sudo apt-get install bwa

.. raw:: html

   <!--
   V: 0.7.5a-r405

   Latest : 0.7.12
   -->

Bowtie2
****************************************************************

`General
documentation <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_

`Instructions <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2>`_

-  Download package
   `here <https://sourceforge.net/projects/bowtie-bio/files/bowtie2/>`_
-  Move package to your personnal bin/
-  Unzip
-  Add to $PATH (see section 1.3)
-  There you go!

::

    cd ~/bin
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip
    unzip bowtie2-2.2.6-linux-x86_64.zip

Peak-calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bPeaks
****************************************************************

Peak-caller developped specifically for yeast, can be useful in order to
process small genomes only.

**TODO**

HOMER
****************************************************************

`Web page <http://homer.salk.edu/>`_

`Install
instructions <http://homer.salk.edu/homer/introduction/install.html>`_

::

    wget "http://homer.salk.edu/homer/configureHomer.pl"
    mkdir $HOME/bin/HOMER
    mv configureHomer.pl $HOME/bin/HOMER
    cd $HOME/bin/HOMERcd $HOME/bin/HOMER
    perl configureHomer.pl -install homer

Add to path (see section 1.3)

::

    PATH="$HOME/bin/HOMER/bin:$PATH"

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
****************************************************************

-  `doc <http://liulab.dfci.harvard.edu/MACS/00README.html>`_
-  `install <http://liulab.dfci.harvard.edu/MACS/INSTALL.html>`_

::

    cd $HOME/bin
    wget "https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz"
    tar -xvzf MACS-1.4.2-1.tar.gz
    cd MACS-1.4.2
    sudo python setup.py install
    macs14 --version

**NB** deb package wouldn't work with python 2.7, asks for python 2.6.

MACS2
****************************************************************

-  `MACS2 web page <https://github.com/taoliu/MACS/>`_

::

    sudo apt-get install python-numpy
    sudo pip install MACS2

.. raw:: html

   <!--
   Marche pas?
   ```
   $ git clone https://github.com/taoliu/MACS.git
   # pip install MACS2
   ...
   ```
   -->

SPP R package (broken)
****************************************************************

::

    install.packages("caTools")
    install.packages("spp")

<!--

::

    source("http://bioconductor.org/biocLite.R")
    biocLite("spp")
    > install.packages("spp")

::

    R CMD INSTALL spp_1.10.tar.gz

...

::

    sudo su
    echo "deb http://www.stats.bris.ac.uk/R/bin/linux/ubuntu precise/" >> /etc/apt/sources.list
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
    apt-get update
    apt-get upgrade

::

    wget "https://cran.r-project.org/src/base/R-3/R-3.2.2.tar.gz"
    tar -xf rm R-3.2.2.tar.gz
    rm R-3.2.2.tar.gz
    cd rm R-3.2.2
    ./configure

doesn't work on VM

not a problem of R version anyway

libboost libraries ? apt-get install libboost-all-dev -->

SWEMBL
****************************************************************

-  `SWEMBL beginner's
   manual <http://www.ebi.ac.uk/~swilder/SWEMBL/beginners.html>`_

**TODO**

Motif discovery, motif analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RSAT suite
****************************************************************

See `doc/install\_protocols
section <https://github.com/rioualen/gene-regulation/blob/master/doc/install_protocols/install_rsat_ubuntu14.04.Rmd>`_.
Beware, this manuel might be deprecated.



Bazar à trier
^^^^^^^^^^^^^^^^^^^^^^^^


Table of Contents

-  `Pre-requisites in case of virtual machine (VM)
   development <#pre-requisites-in-case-of-virtual-machine-vm-development>`_

   -  `VM creation <#vm-creation>`_
   -  `VM customization <#vm-customization>`_

-  `General requirements <#general-requirements>`_

   -  `Generic tools <#generic-tools>`_
   -  `ssh <#ssh>`_
   -  `rsync <#rsync>`_
   -  `git <#git>`_
   -  `zlib <#zlib>`_
   -  `qsub <#qsub>`_
   -  `Create bin/ (opt) <#create-bin-opt>`_
   -  `Edit $PATH <#edit-path>`_

-  `Snakemake workflows basic
   requirements <#snakemake-workflows-basic-requirements>`_

   -  `Python <#python>`_
   -  `Pandas library <#pandas-library>`_
   -  `Package rpy2 <#package-rpy2>`_
   -  `R (to be revised) <#r-to-be-revised>`_
   -  `Snakemake <#snakemake>`_
   -  `Graphviz <#graphviz>`_

-  `NGS analysis software & tools <#ngs-analysis-software--tools>`_

   -  `File management <#file-management>`_
   -  `SRA toolkit <#sra-toolkit>`_
   -  `Samtools <#samtools>`_
   -  `Bedtools <#bedtools>`_
   -  `Quality assessment <#quality-assessment>`_
   -  `FastQC <#fastqc>`_
   -  `Trimming <#trimming>`_
   -  `Sickle <#sickle>`_
   -  `Alignment/mapping <#alignmentmapping>`_
   -  `BWA <#bwa>`_
   -  `Bowtie2 <#bowtie2>`_
   -  `Peak-calling <#peak-calling>`_
   -  `bPeaks <#bpeaks>`_
   -  `HOMER <#homer>`_
   -  `MACS 1.4 <#macs-14>`_
   -  `MACS2 <#macs2>`_
   -  `SPP R package (broken) <#spp-r-package-broken>`_
   -  `SWEMBL <#swembl>`_
   -  `Motif discovery, motif
      analysis <#motif-discovery-motif-analysis>`_
   -  `RSAT suite <#rsat-suite>`_

-  `Workpackage 2.6 - Gene
   regulation <#workpackage-26---gene-regulation>`_

   -  `Cloning the repository <#cloning-the-repository>`_
   -  `Data transfer/download <#data-transferdownload>`_
   -  `Running the pipeline <#running-the-pipeline>`_

-  `VM export / submission <#vm-export--submission>`_ --> <!-- ###
   **TODO**

-  Include map of possible "bricks" of worflows (like general rulegraph)
   with each step's requirement/dependencies
-  Include minimum json file config depending on bricks to be used.

-  Beware of paths

   -  ~/workspace
   -  ~/bin

-  Test all of this in IFB appliance

-  **VBox issues**:

   -  '/etc/init.d/vboxdrv setup' error when restarting
   -  disparition vboxnet0

**/!\\** attention pour les rsync, notamment si user + root... ssh, ssh
agent ?

-  revoir les install via apt-get car pb de version !
-  mettre la procédure spécifique
-  lister les versions de chaque programme pour un wf qui fonctionne

-  dependance mkvtree / rsat

-  see differences between ubuntu and lmde (python libs notamment)
-  check mac ?

-  check version dependencies and add --version to doc -->

