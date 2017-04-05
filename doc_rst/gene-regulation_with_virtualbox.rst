<!--
====

title: "Running the gene-regulation pipelines with VirtualBox" author:
"Claire Rioualen" date: "April 25, 2016" output: html\_document:
fig\_caption: yes highlight: zenburn number\_sections: yes theme:
cerulean toc: yes toc\_depth: 5 pdf\_document: fig\_caption: yes
highlight: zenburn number\_sections: yes toc: yes toc\_depth: 5 --- -->

-  `Creating a virtual machine (VM) <#creating-a-virtual-machine-vm>`__

   -  `Creating a VM under VirtualBox
      software <#creating-a-vm-under-virtualbox-software>`__

      -  `Requirements <#requirements>`__
      -  `Virtual Box configuration <#virtual-box-configuration>`__
      -  `Creation of the virtual
         machine <#creation-of-the-virtual-machine>`__
      -  `VM configuration <#vm-configuration>`__
      -  `Operating system
         installation <#operating-system-installation>`__

-  `Installing programs and
   dependencies <#installing-programs-and-dependencies>`__

   -  `Get the ``gene-regulation``
      repository <#get-the-gene-regulation-repository>`__
   -  `Run makefile to install all required
      dependencies <#run-makefile-to-install-all-required-dependencies>`__

-  `Executing snakemake workflow
   example <#executing-snakemake-workflow-example>`__

   -  `Download source data <#download-source-data>`__
   -  `Execute workflow <#execute-workflow>`__

-  `Visualizing results <#visualizing-results>`__

   -  `FastQC <#fastqc>`__
   -  `IGV <#igv>`__

-  `Export appliance (todo) <#export-appliance-todo>`__
-  `Import appliance (todo) <#import-appliance-todo>`__

    Notes: This protocol was developed on a Unix computer, with the OS
    LMDE and a 64-bit architecture. The virtual machines are developed
    with the Ubuntu 14.04 OS. The protocol is based on gene-regulation
    v2.0

Creating a virtual machine (VM)
===============================

Creating a VM under VirtualBox software
---------------------------------------

Requirements
~~~~~~~~~~~~

**Virtualbox software**

We used VirtualBox 5.0.2, downloadable from https://www.virtualbox.org/
or to be installed manually:

::

    sudo apt-get install virtualbox-5.0

VirtualBox extension pack can be requested (eg. for handling USB2.0, see
'errors' section).

::

    wget http://download.virtualbox.org/virtualbox/5.0.2/Oracle_VM_VirtualBox_Extension_Pack-5.0.2.vbox-extpack

**Ubuntu image**

In this tutorial we used Ubuntu 14.04.4, latest long-term supported
version.

::

    wget http://releases.ubuntu.com/14.04/ubuntu-14.04.4-desktop-amd64.iso

--------------

Virtual Box configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

Before configuring the virtual machine, we need to tell VirtualBox how
it will enable your local virtual machines to interact with their host
(the operating system of the machine on which the VM is running).

1. Open *VirtualBox > File > Preferences...*

2. Open the tab *Network* > *Host-only Networks*

   -  click on the "+" icon
   -  this creates a network vboxnet0. Select this network, click on the
      screw driver icon (*edit host-only network*), and set the
      following options:

   -  *Adapter* tab

      -  IPv4 Address: 192.168.56.1
      -  IPv4 Network Mask: 255.255.255.0
      -  IPv6 Adress: blank
      -  IPv6 Network Mask Length: 0

   -  *DHCP Server* tab

      -  Check *Enable Server*
      -  *Server Address:* 192.168.56.100
      -  *Server Mask:* 255.255.255.0
      -  *Lower Address Bound:* 192.168.56.101
      -  *Upper Address Bound:* 192.168.56.254

.. figure:: ../../img/vbox_network.png
   :alt: 

.. figure:: ../../img/vbox_network_adapter.png
   :alt: 

.. figure:: ../../img/vbox_network_DHCP.png
   :alt: 

Creation of the virtual machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Open VirtualBox

2. Click on the **New** button.

3. Parameters

-  Name and operating system

   -  Name: gene-regulation
   -  Type: Linux
   -  Version: Ubuntu (64 bits)

-  Memory size: 2048 Mb (this can be modified afterwards).

-  Hard drive: *Create a virtual hard drive now*.

-  Hard drive file type: *VDI* (VirtualBox Disk Image).

.. raw:: html

   <!--
   - Hard drive file type: *VMDK* (Virtual Machine Disk). 

       I chose this option because it ensures a wider compatibility with other OS and Virtual Machine management systems. 

          Another potential advantage of VMDK is that it enables to split virtualdisks in files <=2Gb, which is convenient to store them on FAT partitions. 
   -->

-  Storage on physical hard drive

   -  Select *Dynamically allocated*

.. raw:: html

   <!--
       - Activate the option *Split into files less than 2Gb*, which allows to store the VM on FAT partitions for Windows host machines.
   -->

-  File location and size

   -  max size of virtual hard drive: 30GB
   -  click on **Create** button

*Note:* you should adapt the virtual hard drive size to your needs. Be
aware that it's difficult to extend later on, so you should aim larger
than expected. Since the size is dynamically allocated, it won't take up
too much space until you fill it.

At this stage, the VM has been created and needs to be configured before
installing the operating system.

VM configuration
~~~~~~~~~~~~~~~~

In the VirtualBox main window, select the newly created virtual machine,
and click on the **Settings** button.

**General**

For the desktop version of Ubuntu, it is convenient to enable copy-paste
between the guest and the host.

-  Select the tab *Advanced*
-  Set *Shared clipboard* to *Bidirectional*

**Storage**

Click on the **Empty** disc icon in the storage tree. Select the disc
icon on the right and fetch the downloaded ``.iso`` image(see
**Requirements**). Click on *OK*.

**Network**

VirtualBox offers many alternative ways to configure network
communications between the virtual machine, the host machine, and the
external network.

To get more information about network settings:

-  VirtualBox `manual
   page <https://www.virtualbox.org/manual/ch06.html>`__
-  An excellent
   `tutorial <http://christophermaier.name/blog/2010/09/01/host-only-networking-with-virtualbox>`__

We present here one possible way to configure your Virtual machine, but
this should be adapted to the particular security/flexibility
requirements of the network where the maching has to run.

In the VM settings, select tne *Network* tab. VirtualBox enables you to
specify several adapters, each corresponding to one separate network
access (e.g. using an ethernet card + wi-fi connection).

-  click on the tab *Adapter 1*,

   -  check *Enable Network Adapter*
   -  Attached to: *Host-only Adapter*
   -  Name: *vboxnet0* (this network must have been created beforehand,
      see section 1.2.3)

-  click on the tab *Adapter 2*,

   -  check *Enable Network Adapter*
   -  Attached to : *NAT*

-  click on the tab *Adapter 3*,

   -  check *Enable Network Adapter*
   -  Attached to : *Bridged Adapter*
   -  Name: choose an option corresponding to the actual internet
      connection of the host machine (e.g. ethernet cable, Wi-Fi, ...).

**You can now start the VM. **

.. raw:: html

   <!-- This can raise several errors, if so see dedicated section below.  -->

Operating system installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Welcome

   -  check the language settings and click on *Install Ubuntu*.

-  Preparing to install Ubuntu

   -  leave all default parameters and click *Continue*.

-  Installation type

   -  (leave the default) Erase disk and install Ubuntu, click *Install
      Now*.

-  Where are you (automatic)

   -  Paris

-  Keyboard layout

   -  French - French

-  Who are you ?

   -  Your name: gene-regulation
   -  Your computer's name: gene-regulation-virtual
   -  Pick a username: gr
   -  Choose a password: genereg
   -  (Activate the option Log in automatically)

Restart once installation is completed.

Once on the desktop, go to the VM menu: select *Devices* then *Install
Guest Additions CD image*. Run it.

The VirtualBox Guest Additions will provide closer integration between
host and guest and improve the interactive performance of guest systems.
Reboot again to see the new display.

Installing programs and dependencies
====================================

Once in the virtual machine, you can install the required programs from
a terminal.

Get the ``gene-regulation`` repository
--------------------------------------

::

    cd
    wget --no-clobber https://github.com/rioualen/gene-regulation/archive/2.0.tar.gz -P
    tar zvxf 2.0.tar.gz

::

    cd
    git clone https://github.com/rioualen/gene-regulation.git

Run makefile to install all required dependencies
-------------------------------------------------

This may take a while (30mn to 1h) & source the ``.bashrc`` (it's been
updated with the ``$PATH`` for newly installed applications).

::

    cd
    #make -f gene-regulation-2.0/scripts/makefiles/install_tools_and_libs.mk all
    make -f gene-regulation/scripts/makefiles/install_tools_and_libs.mk all
    source ~/.bashrc

Executing snakemake workflow example
====================================

::

    ## Create a base directory for the analysis

    export ANALYSIS_DIR="$HOME/ChIP-seq_SE_GSM20870"
    mkdir $ANALYSIS_DIR

::

    ## Download source data

    mkdir -p ${ANALYSIS_DIR}/data/GSM521934 ${ANALYSIS_DIR}/data/GSM521935
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra -P ${ANALYSIS_DIR}/data/GSM521934
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra -P ${ANALYSIS_DIR}/data/GSM521935

::

    ## Download reference genome & annotations

    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.30.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

::

    ## Execute workflow

    cd ${ANALYSIS_DIR}
    ln -s  $HOME/gene-regulation
    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/ChIP-seq_workflow_SE.py --configfile gene-regulation/examples/ChIP-seq_SE_GSE20870/config.yml

Congratulations! You just executed this wonderful workflow:

.. figure:: ../../img/rule.png
   :alt: 

Visualizing results
===================

FastQC
------

You can visualize the FastQC results using firefox or any other
navigator. Fetch the ``html`` files located in the sample directories.

-  Before trimming:

   ::

       firefox ~/GSE20870-analysis/results/samples/GSM521934/GSM521934_fastqc/GSM521934_fastqc.html
       firefox ~/GSE20870-analysis/results/samples/GSM521935/GSM521935_fastqc/GSM521935_fastqc.html

-  After trimming:

   ::

       firefox ~/GSE20870-analysis/results/samples/GSM521934/GSM521934_sickle-se-q20_fastqc/GSM521934_sickle-se-q20_fastqc.html
       firefox ~/GSE20870-analysis/results/samples/GSM521935/GSM521935_sickle-se-q20_fastqc/GSM521935_sickle-se-q20_fastqc.html

.. figure:: ../../img/vbox_fastqc.png
   :alt: 

IGV
---

You can visualize the peaks by running IGV from the terminal.

::

    igv

-  Click "File" > "Open session..." and chose the file
   ``~/GSE20870-analysis/results/peaks/igv_session.xml``.
-  You may need to adjust the panel sizes.

.. figure:: ../../img/igv.png
   :alt: 

Export appliance (todo)
=======================

The virtual machine created with VirtualBox can be exported and saved as
an appliance.

-  Shut down the VM.
-  In VirtualBox, open *File* -> *Export Appliance ...*

-  Select the VM ``gene-regulation``
-  *Next >*

-  Save as: gene-regulation-[YYMMDD].ova
-  Format: OVF 1.0
-  Write Manifest File: check
-  *Next >*

-  Appliance Settings

   -  Name: gene-regulation-[YYMMDD]
   -  Product: Regulatory Genomics Pipeline
   -  Product-URL: -
   -  Vendor: Claire Rioualen, Jacques van Helden
   -  Version: YYYY-MM-DD
   -  Description: Regulatory Genomics Pipeline using Snakemake,
      installed on an Ubuntu 14.04 Virtual Machine.
   -  License: Free of use for academic users, non-commercial and
      non-military usage.

-  *Export*

The appliance saved can be re-imported later on, on another computer if
needed.

Import appliance (todo)
=======================

In VirtualBox, click menu File > Import appliance > fetch OVA file.

Note: there is apparently a bug with the export of VMs under VirtualBox
5.0. If you get this error when launching the imported file:

    A new node couldn't be inserted because one with the same name
    exists. (VERR\_CFGM\_NODE\_EXISTS).

There is a workaround: go to the imported VM settings, to the USB tab,
and untick "enable USB Controller". You should now be able to start the
VM.
