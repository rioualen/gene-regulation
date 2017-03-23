1. IFB
-----------

Table of Contents (to be removed)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `IFB cloud utilities <#ifb-cloud-utilities>`__

   -  `User account creation &
      configuration <#user-account-creation--configuration>`__

-  `Using an existing appliance <#using-an-existing-appliance>`__

   -  `Virtual disk creation <#virtual-disk-creation>`__
   -  `Creation of an instance <#creation-of-an-instance>`__

-  `Creation of an appliance <#creation-of-an-appliance>`__

   -  `Creation <#creation>`__
   -  `Configuration (optional) <#configuration-optional>`__
   -  `User account <#user-account>`__
   -  `Shell coloring <#shell-coloring>`__
   -  `Data management <#data-management>`__
   -  `Virtual disk <#virtual-disk>`__
   -  `Transfers <#transfers>`__

-  `Then... <#then>`__

   -  `Software installation <#software-installation>`__

IFB cloud utilities
~~~~~~~~~~~~~~~~~~~~

.. figure:: ../../img/ifb-logo.png
   :alt: 

The French Bioinformatics Institute (IFB) cloud provides users with a
number of bioinformatics facilities, under the form of ready-to-use
*appliances*. A cloud appliance is a template or a virtual machine (VM)
built with a bundle of scientific or utility software components that
are already configured. Several appliances are dedicated to special
fields of bioinformatics, such as proteomics, genomics... Some of them
come with an HTML interface, such as Galaxy or RSAT.

The cloud also provides "basic" Ubuntu or CentOS appliances. Provided
you hold a developper account, it allows you to instantiate a virtual
machine, setup your own tools, and register it as a new appliance to be
used again later on and even shared with other cloud users.

The official website is still under development. However, here are a few
useful links:

-  `The IFB <http://www.france-bioinformatique.fr/>`__

-  `IFB cloud <http://www.france-bioinformatique.fr/en/cloud/>`__

-  `Cloud
   usage <http://www.france-bioinformatique.fr/en/core/cloud-usage>`__

-  `Documentation <http://www.france-bioinformatique.fr/en/cloud/doc-du-cloud>`__

User account creation & configuration
**************************************

-  Using the IFB cloud facilities requires to have a user account.
   Register
   `here <https://cloud.france-bioinformatique.fr/accounts/register/>`__.

-  Once your account has been validated, you can
   `login <https://cloud.france-bioinformatique.fr/accounts/login/>`__.

-  In order to be able to access your instances through SSH, you should
   register your SSH public key in your `account
   settings <https://cloud.france-bioinformatique.fr/cloud/profile/>`__,
   through the dashboard.

.. figure:: ../../img/dashboard.png
   :alt: 

Using an existing appliance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Virtual disk creation
***********************

Appliances usually have a limited amount of disk space (up to 10, 20Go).
If the instance to be run necessitates disk space, you have to create a
virtual disk (vDisk) prior to launching it. By default, the capacity of
storage granted to a user is 250Go, which can be divided into as many
vDisks as necessary. When instantiating an appliance, you can chose to
attach one of these vDisks to the virtual machine. You'll be able to
access data on this disk through SSH.

1. Click *New vDisk* button.
2. Enter a size (whole number equating to the amount of Go needed).
3. Name it.

.. figure:: ../../img/create_vDisk.png
   :alt: 

Creation of an instance
**************************

1. Click *New Instance* button.
2. Choose an appliance in the drop-down menu. You may use the filter
   menu in order to look for a specific tool.
3. Name your VM.
4. Choose the amount of CPU and RAM to grant the VM (up to 8 CPU, 32 GB
   RAM).
5. Attach the vDisk.
6. Click *Run*.

.. figure:: ../../img/create_instance.png
   :alt: 

7. After a few seconds, you may refresh the page until the newly created
   instance shows up on the dashboard. Clicking on the ssh mention in
   the *Access* column will give you the commands to access your virtual
   machine.

.. figure:: ../../img/ssh.png
   :alt: 

8. If the appliance has an HTTP interface, a link will also be provided
   in the *Access* column.

Creation of an appliance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creation
*************

Creating your own appliance can be as simple as instantiating an
existing one.

1. Click *New Instance* button.
2. Choose the appliance **Ubuntu 14.04 IFB-10G (2015-10)** or **CentOS
   6.7 IFB-20G (2016-01)**.
3. Name your instance.
4. Check **Create appliance**.
5. Choose the amount of CPU and RAM to grant the VM (up to 8 CPU, 32 GB
   RAM).
6. Attach the vDisk.
7. Click *Run*.

.. figure:: ../../img/create_appliance.png
   :alt: 

8. Refresh the page. Your instance should appear in orange because of
   the creation mode you selected. You can now click on the **ssh**
   column to see the ssh parameters. It should look like this:

.. figure:: ../../img/ubuntu_create.png
   :alt: 

9. Connect to your VM by commandline.

   ::

       ssh -A -p 22 root@192.54.201.111

Configuration (optional)
***************************

User account
^^^^^^^^^^^^^^

Create user account and grant it sudo privileges (followed procedure
`here <https://www.digitalocean.com/community/tutorials/how-to-add-and-delete-users-on-an-ubuntu-14-04-vps>`__).

Shell coloring
^^^^^^^^^^^^^^

::

    nano ~/.bashrc

Fetch following paragraph and uncomment command ``force-color``.

::

    # uncomment for a colored prompt, if the terminal has the capability; turned
    # off by default to not distract the user: the focus in a terminal window
    # should be on the output of commands, not on the prompt
    force_color_prompt=yes

::

    source ~/.bashrc

Data management
**********************

Virtual disk
^^^^^^^^^^^^

By default, if a vDisk has been attached to the VM, it is mounted under
``/root/mydisk``.

Transfers
^^^^^^^^^

You can transfer data from your local computer to the VM using commands
provided under *Access* > ssh:

::

    scp -P 22 ${localfile} root@192.54.201.111:
    sftp -oPort=22 root@192.54.201.111

Another way is to use rsync:

::

    rsync -ruptvl ${localfile} root@192.54.201.177:/root/mydisk/

Then...
~~~~~~~~~~~~~

Software installation
**************************

Once you're connected to the VM through ``ssh``, you can install any
program just the way you would do it locally (see tutorials in `this
directory <https://github.com/rioualen/gene-regulation/tree/master/doc/install_protocols>`__
for instance).

IFB bis et ter

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.1 Requirements
************************

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
****************************

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
****************************

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
*********************************

Open a terminal on your host computer and type in:

::

    ssh -A -p 22 root@192.54.201.124

1.5 Download source data
****************************

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
^^^^^^^^^^^^^^^^

.. raw:: html

   <!--mkdir -p ${ANALYSIS_DIR}/data/GSM521934 ${ANALYSIS_DIR}/data/GSM521935-->

::

    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra -P ${ANALYSIS_DIR}/data/GSM521934
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra -P ${ANALYSIS_DIR}/data/GSM521935

Download reference genome & annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.30.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

You should now have something like this:

.. figure:: ../../img/data_tuto.png
   :alt: 

1.6 Run the workflow
***********************

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

2.1 Install and run the X2Go client on your host computer
************************************************************

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
*************************

The result files should be organized like this:

.. figure:: ../../img/results_orga.png
   :alt: 

FastQC
^^^^^^

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
^^^^^^

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creating a new appliance from scratch is very similar to using one. You
have to satisfy the requirements described in part 1.1.

If you want to manipulate data, you should also create a vDisk following
step 1.2.

Creation of an instance
****************************

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~éé

Once in the virtual machine, you can install the required programs.

Get the ``gene-regulation`` repository
*******************************************

::

    wget https://github.com/rioualen/gene-regulation/archive/2.0.tar.gz
    tar zvxf 2.0.tar.gz

Run makefile to install the dependencies
********************************************

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


2. DOCKER
-------------

<!-- output: html\_document: fig\_caption: yes highlight: zenburn
number\_sections: yes theme: cerulean toc: yes toc\_depth: 5
pdf\_document: fig\_caption: yes highlight: zenburn number\_sections:
yes toc: yes toc\_depth: 5 --- -->

-  `Get started with Docker! <#get-started-with-docker>`__

   -  `Create a Docker account <#create-a-docker-account>`__
   -  `Install Docker on your local
      host <#install-docker-on-your-local-host>`__
   -  `Create shared repositories and download source
      data <#create-shared-repositories-and-download-source-data>`__
   -  `Fetch the Docker image and run it with shared
      folders <#fetch-the-docker-image-and-run-it-with-shared-folders>`__
   -  `Execute the pipeline <#execute-the-pipeline>`__

-  `JVH / Mac <#jvh--mac>`__

    Note: this is a first draft. Tested in Ubuntu 14.04 64bits version.

Get started with Docker!
~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a Docker account
************************

Instructions `here <https://docs.docker.com/linux/step_five/>`__.

Install Docker on your local host
***************************************

Instructions for a linux install can be found
`here <https://docs.docker.com/linux/>`__, along with mac and windows
instructions. A useful script is availalable
`here <https://gist.github.com/bhgraham/ed9f8242dc610b1f38e5>`__ for a
debian install.

You can also install it on Ubuntu 14.04 (64bits) typing the following:

::

    #sudo apt-get update
    sudo apt-get -y install docker.io
    sudo usermod -aG docker <username>

You should now log out and in again from your Ubuntu session. This short
procedure was tested in a virtual machine under VirtualBox (see
corresponding tutorial).

.. raw:: html

   <!--sudo service docker start-->

You can test whether docker works properly:

::

    docker run hello-world

.. figure:: ../../img/docker_hello.png
   :alt: 

NB: it seems qwerty keyboard keeps popping up after docker install.
Switch back to azerty:

::

    setxkbmap fr

<!-- Run the following command:

::

    sudo apt-get --yes install docker

-->

Create shared repositories and download source data
*******************************************************

In order to execute the study case GSE20870, you should enter the
following commands:

::

    export ANALYSIS_DIR=~/GSE20870-analysis
    mkdir $ANALYSIS_DIR
    cd $ANALYSIS_DIR

::

    mkdir data/GSM521934 
    wget -nc ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra -P data/GSM521934

    mkdir data/GSM521935
    wget -nc ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra -P data/GSM521935

Fetch the Docker image and run it with shared folders
********************************************************

::

    docker pull rioualen/gene-regulation:2.0
    docker run -v $ANALYSIS_DIR:~/GSE20870-analysis -it rioualen/gene-regulation:2.0 /bin/bash

You can share as many folders as desired, using this syntax:
``-v /path/on/host/:/path/on/docker/``.

Execute the pipeline
***********************

::

    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/factor_workflow.py --configfile gene-regulation/examples/GSE20870/GSE20870.yml

<!-- # JVH / Mac

Quick tour
**************

On Mac OSX

1. Install docker

::

        https://docs.docker.com/engine/installation/mac/

2. Open the application Docker Quickstart Terminal. This open a new terminal window and launches the docker daemon.


3. Get the gene-regulation docker


docker pull rioualen/gene-regulation:0.3

4. Check the list of docker images available locally


docker images

5. Start the gene-regulation image. The option ``-it`` specifies the interactive mode, which is necessary to be able using this VM


::

    docker run -it rioualen/gene-regulation:0.3 /bin/bash

You are now in a bash session of a gene-regulation docker. In this
session, you are "root" user, i;e. you have all the administration
rights. You can check this easily:

::

    whoami

6. Check the disks available on this docker


::

    df -h

Currently, your docker can only access its local disk, which comes with
the VM. **Beware**: any data stored on this local disk will be lost when
you shut down the gene-regulation docker.

7. Exit and get back to your gene-regulation container


If you exits your shell session, the docker will still be running.

::

    exit

You are now back to the host terminal.

Check the currently active docker containers (processes).

::

    docker ps -a

Note that you can run several containers of the same image. Each active
container has a unique identifier which appears in the first column when
you run ``docker ps`` (e.g. ``faff5298ef95``). You can re-open a running
container with the command

::

    docker attach [CONTAINER_ID]

where ``[CONTAINER_IDR]`` must be replaced by the actual ID of the
running docker container (e.g. ``faff5298ef95``).

8. Shutting down the container


We will now shut down this image, and start a new one which will enable
you to store persistent data.

::

    docker stop [CONTAINER_ID]

9. Starting a docker container with a shared folder.


500 docker pull rioualen/gene-regulation:0.3 501 mkdir -p
~/gene-regulation\_data/GSE20870/GSM521934
~/gene-regulation\_data/GSE20870/GSM521935 502 cd
~/gene-regulation\_data/GSE20870/GSM521934 503 wget
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra
504 cd ~/gene-regulation\_data/GSE20870/GSM521935 505 wget
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra
506 mkdir ~/gene-regulation\_data/results/GSE20870 507 mkdir -p
~/gene-regulation\_data/results/GSE20870 508 docker pull
rioualen/gene-regulation:0.3 509 docker run -v
~/gene-regulation\_data:/data -it rioualen/gene-regulation:0.3 /bin/bash

10. Running the snakemake demo workflow on the docker container


::

    1  ls /data
    2  ls /data/GSE20870/
    3  ls /data/GSE20870/GSM521934/
    4  exit
    5  ls /data
    6  source ~/bin/ngs_bashrc
    7  snakemake -s scripts/snakefiles/workflows/factor_workflow.py -np
    8  history
    9  snakemake -s scripts/snakefiles/workflows/factor_workflow.py -np

10 snakemake -s scripts/snakefiles/workflows/factor\_workflow.py -p

Questions
***********

1. Quand on fait un login dans la vm gene--regulation, on entre dans un
   shell basique (pas bash). Est-il possible de configurer docker pour
   qu'on entre automatiquement en bash ?

Entry point /bin/bash

2. Il faut ajouter le bashrc dans le /etc du docker.

-->



3. VBOX
--------


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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creating a VM under VirtualBox software
********************************************

Requirements
^^^^^^^^^^^^^

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



Virtual Box configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once in the virtual machine, you can install the required programs from
a terminal.

Get the ``gene-regulation`` repository
****************************************

::

    cd
    wget --no-clobber https://github.com/rioualen/gene-regulation/archive/2.0.tar.gz -P
    tar zvxf 2.0.tar.gz

::

    cd
    git clone https://github.com/rioualen/gene-regulation.git

Run makefile to install all required dependencies
**************************************************

This may take a while (30mn to 1h) & source the ``.bashrc`` (it's been
updated with the ``$PATH`` for newly installed applications).

::

    cd
    #make -f gene-regulation-2.0/scripts/makefiles/install_tools_and_libs.mk all
    make -f gene-regulation/scripts/makefiles/install_tools_and_libs.mk all
    source ~/.bashrc

Executing snakemake workflow example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~

FastQC
********

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
******

You can visualize the peaks by running IGV from the terminal.

::

    igv

-  Click "File" > "Open session..." and chose the file
   ``~/GSE20870-analysis/results/peaks/igv_session.xml``.
-  You may need to adjust the panel sizes.

.. figure:: ../../img/igv.png
   :alt: 

Export appliance (todo)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~

In VirtualBox, click menu File > Import appliance > fetch OVA file.

Note: there is apparently a bug with the export of VMs under VirtualBox
5.0. If you get this error when launching the imported file:

    A new node couldn't be inserted because one with the same name
    exists. (VERR\_CFGM\_NODE\_EXISTS).

There is a workaround: go to the imported VM settings, to the USB tab,
and untick "enable USB Controller". You should now be able to start the
VM.


4. CONDA
---------
