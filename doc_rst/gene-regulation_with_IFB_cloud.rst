<!--output: html\_document: fig\_caption: yes highlight: zenburn
number\_sections: yes theme: cerulean toc: yes toc\_depth: 4 toc\_float:
yes pdf\_document: fig\_caption: yes highlight: zenburn
number\_sections: yes toc: yes toc\_depth: 4 word\_document: toc: yes
toc\_depth: 4 --- -->

-  `1. Using the Gene-regulation
   appliance <#1-using-the-gene-regulation-appliance>`__

   -  `1.1 Requirements <#11-requirements>`__
   -  `1.2 Virtual disk creation <#12-virtual-disk-creation>`__
   -  `1.3 Creation of an instance <#13-creation-of-an-instance>`__
   -  `1.4 Connection to the device <#14-connection-to-the-device>`__
   -  `1.5 Download source data <#15-download-source-data>`__
   -  `1.6 Execute workflow <#16-execute-workflow>`__

-  `2. Visualizing results <#2-visualizing-results>`__

   -  `2.1 Install and run the X2Go client on your host
      computer <#21-install-and-run-the-x2go-client-on-your-host-computer>`__
   -  `2.2 Visualize results <#22-visualize-results>`__

      -  `FastQC <#fastqc>`__
      -  `IGV <#igv>`__

-  `3. Create your own Gene-regulation
   appliance <#3-create-your-own-gene-regulation-appliance>`__

   -  `Creation of an instance <#creation-of-an-instance>`__
   -  `Installing programs and
      dependencies <#installing-programs-and-dependencies>`__
   -  `Get the ``gene-regulation``
      repository <#get-the-gene-regulation-repository>`__
   -  `Run makefile to install the
      dependencies <#run-makefile-to-install-the-dependencies>`__

    Notes: This protocol was developed on a Unix computer, with the OS
    LMDE. The virtual machines are developed with the Ubuntu 14.04 OS.
    The protocol is based on gene-regulation v2.0

1. Using the Gene-regulation appliance
======================================

1.1 Requirements
----------------

**User account creation & configuration**

-  Using the IFB cloud facilities requires to have a user account.
   Register
   `here <https://cloud.france-bioinformatique.fr/accounts/register/>`__.

-  Once your account has been validated, you can
   `login <https://cloud.france-bioinformatique.fr/accounts/login/>`__.

-  In order to be able to access your instances through SSH, you should
   register your SSH public key in your `account
   settings <https://cloud.france-bioinformatique.fr/cloud/profile/>`__,
   through the dashboard.

1.2 Virtual disk creation
-------------------------

Appliances usually have a limited amount of disk space (up to 10, 20Go).
If the instance to be run necessitates disk space, you have to create a
virtual disk (vDisk) prior to launching it. By default, the capacity of
storage granted to a user is 250Go, which can be divided into as many
vDisks as necessary. When instantiating an appliance, you can chose to
attach one of these vDIsks to the virtual machine. You'll be able to
access data on this disk through SSH.

1. Click *New vDisk* button.
2. Enter a size (whole number equating to the amount of Go needed).
3. Name it (e.g. ``GSE20870-10Gb``, the ID of the Gene Expression
   Omnibus series that will be stored on the virtual drive).

.. figure:: ../../img/vdisk-x2go.png
   :alt: 

.. raw:: html

   <!--\includegraphics[width=250pt]{img/vdisk-x2go.png}-->

1.3 Creation of an instance
---------------------------

1. Click *New Instance* button.
2. Choose appliance "Gene regulation 2.0" in the drop-down menu.
3. Name your VM.
4. Choose the amount of CPU and RAM to grant the VM (up to 8 CPU, 32 GB
   RAM).
5. Attach the vDisk.
6. Click *Run*.

7. After a few seconds, you may refresh the page until the newly created
   instance shows up on the dashboard. Clicking on the ssh mention in
   the *Access* column will give you the commands to access your virtual
   machine.

.. figure:: ../../img/x2go_ssh.png
   :alt: 

1.4 Connection to the device
----------------------------

Open a terminal on your host computer and type in:

::

    ssh -A -p 22 root@192.54.201.124

1.5 Download source data
------------------------

On the IFB cloud VM, the vDisk is automatically attached and mounted by
default under ``/root/mydisk``, or ``~/mydisk``.

Here we create a folder to store the source data files and download them
.

::

    ANALYSIS_DIR=${HOME}/mydisk/GSE20870-analysis

.. raw:: html

   <!--mkdir -p ${ANALYSIS_DIR}/data -->

.. raw:: html

   <!--mkdir -p ${ANALYSIS_DIR}/genome-->

::

    mkdir -p ${ANALYSIS_DIR}
    cd ${ANALYSIS_DIR}
    ln -s ${HOME}/gene-regulation-2.0 gene-regulation

Download data
~~~~~~~~~~~~~

.. raw:: html

   <!--mkdir -p ${ANALYSIS_DIR}/data/GSM521934 ${ANALYSIS_DIR}/data/GSM521935-->

::

    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra -P ${ANALYSIS_DIR}/data/GSM521934
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra -P ${ANALYSIS_DIR}/data/GSM521935

Download reference genome & annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.30.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

You should now have something like this:

.. figure:: ../../img/data_tuto.png
   :alt: 

1.6 Run the workflow
--------------------

You can use the option ``-n`` to make a dry run.

::

    cd  ${ANALYSIS_DIR}
    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/factor_workflow.py --configfile gene-regulation/examples/GSE20870/GSE20870.yml -n

::

    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/factor_workflow.py --configfile gene-regulation/examples/GSE20870/GSE20870.yml

Using 4CPU & 8Go of RAM, the workflow took about 12mn to complete.

Congratulations! You just executed this wonderful workflow:

.. figure:: ../../img/rule.png
   :alt: 

2. Visualizing results
======================

2.1 Install and run the X2Go client on your host computer
---------------------------------------------------------

The Virtual Machine created on the IFB cloud doesn't have a graphical
interface, but it contains the X2GO software. We're gonna use it to
create a distant desktop to visualize the results from the host machine.

1. Install the x2go client and launch it from your local computer.

::

    sudo apt-get install x2goclient
    x2goclient

.. raw:: html

   <!--2. Copy your ssh key to the authorized keys of the virtual machine. (**à revoir !!**)

   ```
   cat $HOME/.ssh/id_rsa.pub | ssh root@192.54.201.124 "cat >> .ssh/authorized_keys"
   ```
   -->

2. Create a new session using the Mate desktop.

.. figure:: ../../img/x2goclient_session_create.png
   :alt: 

3. The session now appears on the right panel. Just click it to lauch
   it!

.. figure:: ../../img/x2go_launch_session.png
   :alt: 

4. You should be now on the virtual desktop!

.. figure:: ../../img/mate_term.png
   :alt: 

Note: you may need to change your keyboard settings

-  Go to **System** > **Preferences** > **Keybords**
-  Click on tab **Layouts**
-  Add and/or remove desired keyboards

2.2 Visualize results
---------------------

The result files should be organized like this:

.. figure:: ../../img/results_orga.png
   :alt: 

FastQC
~~~~~~

You can visualize the FastQC results using firefox or any other
navigator. Fetch the ``html`` files located in the sample directories.

-  Before trimming:

   ::

       firefox /root/mydisk/GSE20870-analysis/results/samples/GSM521934/GSM521934_fastqc/GSM521934_fastqc.html
       firefox /root/mydisk/GSE20870-analysis/results/samples/GSM521935/GSM521935_fastqc/GSM521935_fastqc.html

-  After trimming:

   ::

       firefox /root/mydisk/GSE20870-analysis/results/samples/GSM521934/GSM521934_sickle-se-q20_fastqc/GSM521934_sickle-se-q20_fastqc.html
       firefox /root/mydisk/GSE20870-analysis/results/samples/GSM521935/GSM521935_sickle-se-q20_fastqc/GSM521935_sickle-se-q20_fastqc.html

.. figure:: ../../img/x2go_fastqc.png
   :alt: 

IGV
~~~

You can visualize the peaks by running IGV from the terminal.

.. raw:: html

   <!--You may need to source the `~/.bashrc` first in order to update the `$PATH`. 
   ```
   source ~/.bashrc
   -->

::

    igv

-  Click "File" > "Open session..." and chose the file
   ``/root/mydisk/GSE20870-analysis/reports/peaks/igv_session.xml``.
-  You may need to adjust the panel sizes.

.. figure:: ../../img/igv.png
   :alt: 

3. Create your own Gene-regulation appliance
============================================

Creating a new appliance from scratch is very similar to using one. You
have to satisfy the requirements described in part 1.1.

If you want to manipulate data, you should also create a vDisk following
step 1.2.

Creation of an instance
-----------------------

When creating a new instance choose a 10 Go Ubuntu appliance and check
the "Create appliance" option:

1. Click *New Instance* button.
2. **Choose appliance "Ubuntu 14.04 IFB-X2GO-10GB" in the drop-down
   menu.**
3. Name your VM.
4. Choose the amount of CPU and RAM to grant the VM (up to 8 CPU, 32 GB
   RAM).
5. **Check the box *Create appliance*.**
6. Attach the vDisk.
7. Click *Run*.

.. figure:: ../../img/create_appliance.png
   :alt: 

The new instance should appear in orange bold fonts in the dashboard.

.. figure:: ../../img/ubuntu_create.png
   :alt: 

You can connect to the instance through ``ssh`` as shown in part 1.4.

Installing programs and dependencies
====================================

Once in the virtual machine, you can install the required programs.

Get the ``gene-regulation`` repository
--------------------------------------

::

    wget https://github.com/rioualen/gene-regulation/archive/2.0.tar.gz
    tar zvxf 2.0.tar.gz

Run makefile to install the dependencies
----------------------------------------

This may take a while (up to 30mn-1h) & source the ``.bashrc`` in order
to update the ``$PATH`` accordingly.

::

    make -f gene-regulation-2.0/scripts/makefiles/install_tools_and_libs.mk all
    source ~/.bashrc

If you want to install the x2go server on the VM for visualization
purposes:

::

    make -f gene-regulation-2.0/scripts/makefiles/install_tools_and_libs.mk desktop_and_x2go

You should now be able to execute the example workflow by following
steps 1.5 and 1.6.

In order for your appliance to remain persistant and be available to
other users on the IFB cloud, you should contact an admin: @?
