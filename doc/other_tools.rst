Other manuals (*to be updated*)
================================================================


RSAT installation manual
----------------------------------------------------------------

NB: OLD, to be updated / hopefully replaced with apt-get 


Unless stated otherwise, the following commands will be executed as
root.

!!! Experimented users only !!!

First steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set the locales manually if needed.


::

    export LANGUAGE=en_US.UTF-8
    export LANG=en_US.UTF-8
    export LC_ALL=en_US.UTF-8
    locale-gen en_US.UTF-8
    sudo dpkg-reconfigure locales

Check date & time, and adjust your timezone if needed.


::

    date

::

    dpkg-reconfigure tzdata

Create RSAT directory. It must be readable by all users (in particular by the apache user).


::

    mkdir -p /packages
    cd /packages

Check the installation device. This should give sda1 or vda1?


::

    DEVICE=`df -h / | grep '/dev'| awk '{print $1}' | perl -pe 's/\/dev\///'`
    echo ${DEVICE}

Packages installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Update apt-get.


::

    apt-get update
    apt-get --quiet --assume-yes upgrade

Packages to be installed.


::

    PACKAGES="
    ssh
    git
    cvs
    wget
    zip
    unzip
    finger
    screen
    make
    g++
    apache2
    php5
    libapache2-mod-php5
    libgd-tools
    libgd-gd2-perl
    ghostscript
    gnuplot
    graphviz
    mysql-client
    default-jre
    python
    python-pip
    python-setuptools 
    python-numpy
    python-scipy
    python-matplotlib
    python-suds
    python3
    python3-pip
    python3-setuptools 
    python3-numpy
    python3-scipy
    python3-matplotlib
    r-base
    emacs
    x11-apps
    firefox
    eog
    ntp
    curl
    libcurl4-openssl-dev
    "

Perl modules.


::

    PACKAGES_PERL="perl-doc
    pmtools
    libyaml-perl
    libemail-simple-perl
    libemail-sender-perl
    libemail-simple-creator-perl
    libpostscript-simple-perl
    libstatistics-distributions-perl
    libio-all-perl
    libobject-insideout-perl
    libobject-insideout-perl
    libsoap-lite-perl
    libsoap-wsdl-perl
    libxml-perl
    libxml-simple-perl
    libxml-compile-cache-perl
    libdbi-perl
    liblockfile-simple-perl
    libobject-insideout-perl
    libgd-perl
    libdbd-mysql-perl
    libjson-perl
    libbio-perl-perl
    libdigest-md5-file-perl
    libnet-address-ip-local-perl
    "

Install the apt-get libraries.


::

    echo "Packages to be installed with apt-get --quiet --assume-yes"
    echo "${PACKAGES}"
    echo "Perl module packages to be installed with apt-get --quiet --assume-yes"
    echo "${PACKAGES_PERL}"
    for LIB in ${PACKAGES} ${PACKAGES_PERL}; \
    do \
       echo "`date '+%Y/%m/%d %H:%M:%S'`  installing apt-get library ${LIB}" ; \
       sudo apt-get install --quiet --assume-yes ${LIB} ; \
    done

Package to be installed in an interactive mode.


::

    apt-get install --quiet --assume-yes console-data

-  Options:

   -  Select keymap from arch list
   -  <Don't touch keymap> (default)
   -  Keep kernel keymap
   -  Select keymap from full list

Specific treatment for some Python libraries.


::

    sudo apt-get --quiet --assume-yes build-dep python-numpy python-scipy

To free space, remove apt-get packages that are no longer required. /?\\


::

    apt-get --quiet --assume-yes  autoremove
    apt-get --quiet --assume-yes  clean

Python libraries installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install soappy
    pip install fisher
    pip install httplib2

Apache Web server configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**/!\\** Manual interventions needed here.

Activate CGI module.


::

    nano /etc/apache2/sites-available/000-default.conf

Uncomment the following line:
``Include conf-available/serve-cgi-bin.conf``.

To avoid puzzling warning at apache start, set ServerName globally.


::

    nano /etc/apache2/apache2.conf

Add the following line at the end of the file: ``ServerName localhost``.

Add CGI script.


::

    nano /etc/apache2/mods-available/mime.conf

Uncomment the line ``AddHandler cgi-script .cgi``.

Optional: associate a plain/text mime type to extensions for some
classical bioinformatics files. ``AddType text/plain .fasta``
``AddType text/plain .bed``.

Adapt the PHP parameters.


::

    nano /etc/php5/apache2/php.ini

Modify the following parameters: ``post_max_size = 100M`` and
``upload_max_filesize=100M``.

Activate cgi scripts. Found `here <http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/>`__.


::

    chmod 755 /usr/lib/cgi-bin
    chown root.root /usr/lib/cgi-bin
    a2enmod cgi
    service apache2 restart

You can check whether apache server was successfully configured and
started by opening a web connection to ``http://{IP}``.

RSAT distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**/!\\ Note:** The git distribution requires an account at the ENS git
server, which is currently only possible for RSAT developing team. In
the near future, we may use git also for the end-user distribution. For
users who don't have an account on the RSAT git server, the code can be
downloaded as a tar archive from the Web site.

Create RSAT directory.


::

    mkdir -p /packages/rsat
    cd /packages
    export RSAT=/packages/rsat

Git repository cloning.


::

    git clone git@depot.biologie.ens.fr:rsat
    git config --global user.mail claire.rioualen@inserm.fr
    git config --global user.name "reg-genomics VM user"

\*\* OR \*\*

Archive download.


::

    export RSAT_DISTRIB=rsat_2016-11-06.tar.gz
    export RSAT_DISTRIB_URL=http://pedagogix-tagc.univ-mrs.fr/download_rsat/${RSAT_DISTRIB}

::

    sudo wget ${RSAT_DISTRIB_URL}
    sudo tar -xpzf ${RSAT_DISTRIB}
    sudo rm -f ${RSAT_DISTRIB}
    cd ~; ln -fs /packages/rsat rsat

RSAT configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the configuration script, to specify the environment variables.


::

    cd $RSAT
    sudo perl perl-scripts/configure_rsat.pl

Which options to specify?

.. raw:: html

   <!---
   Absolute path to the RSAT package ? [/packages/rsat]
   Ready to update config file /packages/rsat/RSAT_config.props [y/n] (y):
   rsat_site [your_server_name] : 
   rsat_server_admin [your.mail@your.mail.server] :
   RSAT [/packages/rsat] :
   group_specificity [None] : 
   phylo_tools [0] : 
   variations_tools [0] : 
   ucsc_tools [0] :
   ensembl_tools [0] :
   RSAT_BIN [/packages/rsat/bin] :
   rsat_tmp [/packages/rsat/public_html/tmp] :
   mail_supported [no] : 
   smtp [] : 
   smtp_sender [] : 
   rsat_www [auto] : 
   rsat_echo [0] : 
   start_time [0] : 
   exec_time [0] : 
   rsat_ws [http://localhost/rsat/] : 
   rsat_img_format [png] : 
   QUEUE_MANAGER [batch] : 
   CLUSTER_QUEUE [rsat] : 
   BATCH_MAIL [a] : 
   CLUSTER_SHELL [/bin/bash] : 
   QSUB_OPTIONS [] : 
   REFSEQ_DIR [/packages/rsat/downloads/ftp.ncbi.nih.gov/genomes/refseq] : 
   ensembl_host [ensembldb.ensembl.org] : 
   ensembl_rsync [rsync://ftp.ensembl.org/ensembl/pub] : 
   ensembl_version [79] : 
   ensemblgenomes_version [26] : 
   ensembl_version_safe [70] : 
   ensembl [/packages/rsat/lib/ensemblgenomes-26-79/ensembl/modules] : 
   compara [/packages/rsat/lib/ensemblgenomes-26-79/ensembl-compara/modules] : 
   variation [/packages/rsat/lib/ensemblgenomes-26-79/ensembl-variation/modules] : 
   neat_supported [1] : 
   neat_www_root [http://wwwsup.scmbb.ulb.ac.be/rsat/] : 
   neat_ws [http://wwwsup.scmbb.ulb.ac.be/rsat/web_services/RSATWS.wsdl] : 
   neat_ws_tmp [http://wwwsup.scmbb.ulb.ac.be/rsat/tmp/] : 
   neat_java_ws [http://wwwsup.scmbb.ulb.ac.be/be.ac.ulb.bigre.graphtools.server/wsdl/GraphAlgorithms.wsdl] : 
   neat_java_host [http://wwwsup.scmbb.ulb.ac.be/rsat/] : 
   tomcat_port [] : 
   REA_ROOT [/packages/rsat/contrib/REA] : 
   KWALKS_ROOT [/packages/rsat/contrib/kwalks/bin] : 
   LOGO_PROGRAM [seqlogo] : 
   Ready to update config file /packages/rsat/RSAT_config.mk [y/n] (y): 
   RSAT_SITE [your_server_name] : 
   RSAT_SERVER_ADMIN [your.mail@your.mail.server] : 
   OS [linux] : 
   ARCHITECTURE [x64] : 
   PACKAGE_MANAGER [apt-get] : 
   UCSC_OS [linux.x86_64] : 
   SRC_DIR [${RSAT}/app_sources] : 
   SUDO [] : 
   RSAT_BIN [/packages/rsat/bin] : 
   RSAT_WS [http://localhost/rsat/] : 
   QUEUE_MANAGER [batch] : 
   CLUSTER_QUEUE [rsat] : 
   ENSEMBL_RELEASE [79] : 
   ENSEMBLGENOMES_BRANCH [26] : 
   Ready to update config file /packages/rsat/RSAT_config.bashrc [y/n] (y): 
   Ready to update config file /packages/rsat/RSAT_config.conf [y/n] (y): 
   -->

Load the (updated) RSAT environment variables.


::

    source RSAT_config.bashrc

Check that the RSAT environment variable has been properly configured.


::

    echo ${RSAT}

Initialise RSAT folders


::

    make -f makefiles/init_rsat.mk init

Perl modules for RSAT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cpan

::

    cpan> install YAML
    cpan> install CPAN 
    cpan> reload cpan
    cpan> quit

Get the list of Perl modules to be installed.


::

    make -f makefiles/install_rsat.mk  perl_modules_list
    make -f makefiles/install_rsat.mk perl_modules_check
    more check_perl_modules_eval.txt
    grep Fail  check_perl_modules_eval.txt
    grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;'
    MISSING_PERL_MODULES=`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | cut -f 2 | xargs`
    echo "Missing Perl modules:     ${MISSING_PERL_MODULES}"

Install the missing Perl modules.


::

    make -f makefiles/install_rsat.mk perl_modules_install PERL_MODULES="${MISSING_PERL_MODULES}"

Check once more if all required Perl modules have been correctly installed.


::

    make -f makefiles/install_rsat.mk perl_modules_check
    more check_perl_modules_eval.txt

Note: Object::InsideOut always displays "Fail", whereas it is OK during
installation.

Configure RSAT web server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cd ${RSAT}
    sudo rsync -ruptvl RSAT_config.conf /etc/apache2/sites-enabled/rsat.conf
    apache2ctl restart

RSAT Web server URL


::

    echo $RSAT_WWW

If the value is "auto", get the URL as follows:


::

    export IP=`ifconfig eth0 | awk '/inet /{print $2}' | cut -f2 -d':'`
    echo ${IP}
    export RSAT_WWW=http://${IP}/rsat/
    echo $RSAT_WWW

Other
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

compile RSAT programs written in C


::

    make -f makefiles/init_rsat.mk compile_all
    export INSTALL_ROOT_DIR=/packages/

Install some third-party programs required by some RSAT scripts.


::

    make -f makefiles/install_software.mk install_ext_apps

Mkvtree licence / Vmatch


Get a licence `here <http://www.vmatch.de/>`__

Alternately, you can copy-paste from another RSAT device...

::

    rsync -ruptvl /packages/rsat/bin/vmatch.lic root@<IP>:/packages/rsat/bin/

Data management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    export RSAT_DATA_DIR=/root/mydisk/rsat_data
    cd ${RSAT}/public_html
    mv data/* ${RSAT_DATA_DIR}/
    mv data/.htaccess ${RSAT_DATA_DIR}/
    rmdir data
    ln -s ${RSAT_DATA_DIR} data
    cd $RSAT

Install model organisms, required for some of the Web tools.


::

    download-organism -v 1 -org Saccharomyces_cerevisiae -org Escherichia_coli_K_12_substr__MG1655_uid57779
    download-organism -v 1 -org Drosophila_melanogaster

Get the list of organisms supported on your computer.


::

    supported-organisms

Install selected R librairies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Packages required for some RSAT scripts.


::

    cd $RSAT; make -f makefiles/install_rsat.mk install_r_packages

::

    cd $RSAT; make -f makefiles/install_rsat.mk update ## install R packages + compile the C programs

NB: second only if git repo

Testing RSAT & external programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test a simple Perl script that does not require for organisms to be installed.(OK)


::

    which random-seq
    random-seq -l 100

Test a simple python script that does not require organisms to be installed.(OK)


::

    random-motif -l 10 -c 0.90

Test vmatch


::

    random-seq -l 100 | purge-sequence

seqlogo


::

    which seqlogo
    seqlogo

weblogo 3


::

    which weblogo
    weblogo --help

ghostscript


::

    which gs
    gs --version

Check that the model genomes have been correctly installed


::

    # Retrieve all the start codons and count oligonucleotide frequencies (most should be ATG).
    retrieve-seq -org Saccharomyces_cerevisiae -all -from 0 -to +2 | oligo-analysis -l 3 -1str -return occ,freq -sort

Configure the SOAP/WSDL Web services
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the URL of the web services (RSAT\_WS). By default, the server
addresses the WS requests to itself (http://localhost/rsat) because web
services are used for multi-tierd architecture of some Web tools
(retrieve-ensembl-seq, NeAT).

::

    cd $RSAT
    #echo $RSAT_WS

Get the current IP address


::

    export IP=`/sbin/ifconfig eth0 | awk '/inet /{print $2}' | cut -f2 -d':'`
    echo ${IP}
    export  RSAT_WS=http://${IP}/rsat/

Initialize the Web services stub


::

    make -f makefiles/init_rsat.mk ws_init RSAT_WS=${RSAT_WS}

After this, re-generate the web services stubb, with the following command


::

    make -f makefiles/init_rsat.mk ws_stub RSAT_WS=${RSAT_WS}

Test the local web services OK


::

    make -f makefiles/init_rsat.mk ws_stub_test

Test RSAT Web services (local and remote) without using the SOAP/WSDL stubb (direct parsing of the remote WSDL file)


::

    make -f makefiles/init_rsat.mk ws_nostub_test

Test the program supported-organisms-server, which relies on Web services without stub


::

    supported-organisms-server -url ${RSAT_WS} | wc
    supported-organisms-server -url http://localhost/rsat/ | wc
    supported-organisms-server -url http://rsat-tagc.univ-mrs.fr/ | wc

Tests on the Web site


Run the demo of the following tools (**to redo**)

-  retrieve-seq to check the access to local genomes (at least
   Saccharomyces cerevisiae)
-  feature-map to check the GD library
-  retrieve-ensembl-seq to check the interface to Ensembl
-  fetch-sequences to check the interface to UCSC
-  some NeAT tools (they rely on web services)
-  peak-motifs because it mobilises half of the RSAT tools -> a good
   control for the overall installation.
-  footprint-discovery to check the tools depending on homology tables
   (blast tables).

Install the cluster management system (torque, qsub, ...)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the number of core (processors)


::

    grep ^processor /proc/cpuinfo

Check RAM


::

    grep MemTotal /proc/meminfo

Install Sun Grid Engine (SGE) job scheduler


Beware, before installing the grid engine we need to modify manually the file ``/etc/hosts``


::

    nano /etc/hosts

Initial config (problematic)

::

    127.0.0.1       localhost       rsat-vm-2015-02
    127.0.1.1      rsat-vm-2015-02

Config to obtain:

::

    127.0.0.1       localhost       rsat-vm-2015-02
    127.0.1.1      rsat-vm-2015-02

**/?\\**

::

    apt-get install --quiet --assume-yes gridengine-client
    apt-get install --quiet --assume-yes gridengine-exec
    apt-get install --quiet --assume-yes gridengine-master
    apt-get install --quiet --assume-yes gridengine-qmon 

::

    qconf -aq default  ## aggregate a new queue called "default"
    qconf -mq default  ## modify the queue "default"
    qconf -as localhost ## aggregate the localhost tho the list of submitters

Set the following values: ``hostlist              localhost``

Take all default parameters BUT for the SGE master parameter, type
``localhost`` (it must be the hostname)

Test that jobs can be sent to the job scheduler.

OPTIONAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install some software tools for NGS analysis.


::

    cd ${RSAT}
    make -f makefiles/install_software.mk install_meme

Ganglia: tool to monitor a cluster (or single machine)


`Link. <https://www.digitalocean.com/community/tutorials/introduction-to-ganglia-on-ubuntu-14-04>`__

::

    sudo apt-get install -y ganglia-monitor rrdtool gmetad ganglia-webfrontend
    sudo cp /etc/ganglia-webfrontend/apache.conf /etc/apache2/sites-enabled/ganglia.conf
    sudo apachectl restart



Galaxy server setup 
----------------------------------------------------------------

Downloading Galaxy code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We followed the instructions from the Galaxy Web site:

-  https://wiki.galaxyproject.org/Admin/GetGalaxy

\`\`\`{r eval=FALSE} ## get a git clone of galaxy git clone
https://github.com/galaxyproject/galaxy/ cd galaxy ## Go th the galaxy
directory

Check out the master branch, recommended for production server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| git checkout -b master origin/master
| git pull ## Just in case, we are already up-to-date \`\`\`

Configure the Galaxy server (and get python modules if required)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first edit the config file to chooe a specific port for Galaxy

``{r eval=FALSE} cp config/galaxy.ini.sample config/galaxy.ini``

We then edit this file by setting the port to 8082, because our 8080 is
already used for other purposes.

We performed the following modifications.

admin\_users=admin1@address.fr,admin2@univbazar.fr,admin3@gmail.com port
= 8082 # The port on which to listen. host = 0.0.0.0 ## To enable access
over the network allow\_user\_deletion = True

Configuring the Apache server on RSAT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Activate the Apache module rewrite.load

``{r eval=FALSE} ln -s /etc/apache2/mods-available/rewrite.load  /etc/apache2/mods-enabled/rewrite.load``

Create a file /etc/apache2/sites-enabled/galaxy.conf with the following
content

::

    <VirtualHost *:80>
    ServerAdmin webmaster@localhost
    ServerSignature Off

    # Config pour galaxy ands http://mydomain.com/galaxy
    RewriteEngine on
    RewriteRule ^/galaxy$ /galaxy/ [R]
    RewriteRule ^/galaxy/static/style/(.*) /home/galaxy/galaxy/static/june_2007_style/blue/$1 [L]
    RewriteRule ^/galaxy/static/scripts/(.*) /home/galaxy/galaxy/static/scripts/packed/$1 [L]
    RewriteRule ^/galaxy/static/(.*) /home/galaxy/galaxy/static/$1 [L]
    RewriteRule ^/galaxy/favicon.ico /home/galaxy/galaxy/static/favicon.ico [L]
    RewriteRule ^/galaxy/robots.txt /home/galaxy/galaxy/static/robots.txt [L]
    RewriteRule ^/galaxy(.*) http://localhost:8082$1 [P]
    #RewriteRule ^/galaxy(.*) http://192.168.1.6:8082$1 [P]
    </VirtualHost>

Restart the Apache server.
``{r eval=FALSE} sudo service apache2 restart``

Starting the galaxy server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``{r eval=FALSE} sh run.sh``

On our internal network, the server becomes available at the address:

http://192.168.1.6:8082

Registrating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  open a connection to the Galaxy server
-  In the Galaxy menu, run the command **User -> Register**. Enter the
   same email address as you declared as admin users.

Install Galaxy modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
