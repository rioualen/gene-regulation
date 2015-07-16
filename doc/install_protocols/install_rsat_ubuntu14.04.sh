################################################################
## Instructions used to install a Virtual Machine on the IFB cloud
## (Institut Francais de Bioinformatique), and on a VirtualBox VM.

## THIS IS NOT REALLY A bahsrc FILE, IT IS A SUCCESSION OF
## INSTRUCTIONS AND COMMENTS, THAT SHOULD BE DONE MANUALLY. I (JvH)
## SHOULD IMPROVE THIS WHEN I CAN.

## TO RUN Virtualbox on Mac OSX
#    Go into your Ubuntu Settings
#    Navigate to Keyboard → Keyboard Layout Settings
#    add English (Macintosh)


## For some OS, I must set the locales manually
export LANGUAGE=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
locale-gen en_US.UTF-8
sudo dpkg-reconfigure locales ## as root

################################################################
## Must be executed as root. If you are non-root but sudoer user, you
## can become it withn "sudo bash"

################################################################
## Before anything else, check that the date, time and time zone are
## correctly specified
date

## If not, set up the time zone, date and time with this command
## (source: https://help.ubuntu.com/community/UbuntuTime).
dpkg-reconfigure tzdata

## Configuration for the installation
export INSTALLER=apt-get
export INSTALLER_OPT="--quiet --assume-yes"
## alternative: INSTALLER=aptitude
# export INSTALL_ROOT_DIR=/bio                        ## removed
# export RSAT_HOME=/bio/rsat
# export RSAT_DISTRIB=rsat_2014-08-22.tar.gz
# export RSAT_DISTRIB_URL=http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/${RSAT_DISTRIB}
# 

## We need to update apt-get, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu

## Create a separate directory for RSAT, which must be readable by all
## users (in particular by the apache user)
mkdir -p /bio ## as root
cd /bio
mkdir -p /bio/install_logs ## as root
chmod 777 /bio/install_logs ## as root
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_start.txt


## Check the installation device 
DEVICE=`df -h / | grep '/dev'| awk '{print $1}' | perl -pe 's/\/dev\///'`
echo ${DEVICE}
## This should give sda1

## We can then check the increase of disk usage during the different
## steps of the installation
grep ${DEVICE} /bio/install_logs/df_*.txt

## Install aptitude, more efficient than apt-get to treat dependencies
## when installing and uninstalling packages.
## TO SAVE SPACE, I SUPPRESS aptitude
## apt-get install aptitude
apt-get update ## as root
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_apt-get_updated.txt
${INSTALLER} ${INSTALLER_OPT} upgrade ## as root
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${INSTALLER}_upgraded.txt
grep ${DEVICE} /bio/install_logs/df_*.txt

## Packages to be checked: to I really need this ?
PACKAGES_OPT="
ess
yum
php-elisp
libgdbm-dev
libgd2-xpm-dev
libxml2-dev
links
gfortran
libmysqlclient-dev
texlive-latex-base
python-virtualenv
ipython
ipython-notebook
libssl-dev
libreadline-gplv2-dev:i386
lib64readline-gplv2-dev:i386
libreadline-gplv2-dev
libx11-dev
libxt-dev
libxml2-dev
tcl8.5-dev
tk8.5-dev
libxss-dev
libpng12-dev
libjpeg62-dev
libcairo2-dev
lib32z1
lib32ncurses5
lib32bz2-1.0
libc6-dev
build-essential
python-dev
python3-dev
libnet-ssleay-perl
libcrypt-ssleay-perl
exfat-fuse
exfat-utils 
at
"

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

################################################################
## apt-get packages to install Perl modules (not properly speaking
## necessary, could be done with cpan, but ensure consistency with
## ubuntu OS)

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


## We did not find apt-get packages for some required Perl
## libraries. These will have to be installed with cpan.
PACKAGES_PERL_MISSING="
libalgorithm-cluster-perl
digest-md5-file-perl
liblockfile-simple
libutil-properties-perl
librest-client-perl
libxml-compile-soap11-perl
libxml-compile-wsdl11-perl
libxml-compile-transport-soaphttp-perl
libbio-das-perl        
"

## Install the apt-get libraries ## as root
echo "Packages to be installed with ${INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES}"
echo "Perl module packages to be installed with ${INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES_PERL}"
for LIB in ${PACKAGES} ${PACKAGES_PERL}; \
do \
   echo "`date '+%Y/%m/%d %H:%M:%S'`  installing apt-get library ${LIB}" ; \
   sudo ${INSTALLER} install ${INSTALLER_OPT} ${LIB} > /bio/install_logs/${INSTALLER}_install_${LIB}.txt ; \
   df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${LIB}_installed.txt ; \
done
echo "Log files are in folder /bio/install_logs"
grep ${DEVICE} /bio/install_logs/df_*.txt


## Warning message : 
# apache2_invoke: Enable module perl
# AH00558: apache2: Could not reliably determine the server's fully qualified domain name, using 127.0.1.1. Set the 'ServerName' directive globally to suppress this message

## This package has to be installed in an interactive mode (dialog
## box)
${INSTALLER} install ${INSTALLER_OPT} console-data ## as root

options:
  - Select keymap from arch list
  - <Don't touch keymap> (default)
  - Keep kernel keymap
  - Select keymap from full list

################################################################
## Specific treatment for some Python libraries
##
## A fix for a problem to install scipy with pip: use ${INSTALLER} build-dep 
## taken from here: http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
## Note that these dependencies cost 400Mb ! To be checked
sudo ${INSTALLER} ${INSTALLER_OPT} build-dep python-numpy python-scipy
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_numpy-scipy_dependencies_installed.txt
grep ${DEVICE} /bio/install_logs/df_*.txt

################################################################
## To free space, remove apt-get packages that are no longer required.a
grep ${DEVICE} /bio/install_logs/df_*.txt
${INSTALLER} ${INSTALLER_OPT}  autoremove ## as root
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_autoremoved.txt
${INSTALLER} ${INSTALLER_OPT}  clean ## as root
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned.txt
## This really helps: it saves several hundreds Mb (39% -> 34%?)
grep ${DEVICE} /bio/install_logs/df_*.txt

## DONE: installation of Ubuntu packages
################################################################


################################################################
## Activate the Apache Web server
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!! SOME MANUAL INTERVENTION IS REQUIRED HERE  !!!!!!!!! ## as root
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sudo vi /etc/apache2/sites-available/000-default.conf
## Uncomment the following line:
# Include conf-available/serve-cgi-bin.conf


## To avoid puzzling warning at apache start, set ServerName globally.
sudo nano /etc/apache2/apache2.conf
## Add the following line at the end of the file (or somewhere else)
##     ServerName localhost
                                           

sudo nano /etc/apache2/mods-available/mime.conf
## In the file /etc/apache2/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi
##
## Optional : also associate a plain/text mime type to extensions for
## some classical bioinformatics files.
##   AddType text/plain .fasta
##   AddType text/plain .bed

## Adapt the PHP parameters
sudo nano /etc/php5/apache2/php.ini
## Modify the following parameters
##      post_max_size = 100M
##      upload_max_filesize=100M


## The following lines are required to activate cgi scripts.  Found at
## http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/
sudo chmod 755 /usr/lib/cgi-bin
sudo chown root.root /usr/lib/cgi-bin
sudo a2enmod cgi ## this is apparently required to enable cgi
service apache2 restart

## DONE: apache server configured and started
## You can check it by opening a Web connection to 
## http://[IP]
################################################################

################################################################
## Install some python libraries with pip
##
## Note: at this stage, numpy, scipy and matplotlib have already been
## installed with apt-get under Ubuntu. For other OS, they should be
## added to the pip installation
sudo pip install soappy
sudo pip install fisher
sudo pip install httplib2
## pip install pygraphviz ## OSError: Error locating graphviz.

## optional: an utility to measure internet bandwidth
sudo pip install speedtest-cli

#${INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-suds

## Failures: no distributions at all found
# pip3 install wsdl
# pip3 install wstools
sudo pip3 install fisher
sudo pip3 install snakemake
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## PROBLEM: soappy seems to be discontnued for python3 !
#      pip3 install soappy
## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages ### skipping
pip3 install suds-jurko
pip3 install pysimplesoap

## Check disk usage
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
grep ${DEVICE} /bio/install_logs/df_*.txt

# ################################################################ skipping ########
# ## TO BE CHECKED: TO WE STILL NEED TO DO ALL THE TRICKY STUFF BELOW ?
# ## NOT SURE: on 2014/08/15 Jv installed a VM on IFB cloud, and the
# ## SOAP/WSDL Web services seem to work without this.
#
# ## The installation of SOAP:WSDL under cpan is particularly tricky. 
# ## In Ubuntu, there is a way to install it with ${INSTALLER}. 
# ## http://www.installion.co.uk/ubuntu/trusty/universe/l/libsoap-wsdl-perl/fr/install.html
# emacs -nw /etc/apt/sources.list
#
# ## Ensure that the following line is set to "universe"
# deb http://us.archive.ubuntu.com/ubuntu trusty main universe
# ## You can now quit emacs
#
# apt-get update
#
# apt-get --quiet --assume-yes install libmodule-build-perl
# apt-get --quiet --assume-yes install libsoap-wsdl-perl
#
# ## Note: this is still not sufficient to get SOAP::WSDL to run the two
# ## following targets
# ##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stubb
# ##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stubb_test
#
# ## We first need to fix some problem with CPAN : on Ubuntu, I cannot
# ## install the SOAP::WSDL module, which is required for several
# ## functionalities. More precisely, after fiddling around for a few
# ## hours, the server is able to answer to web services requests, but I
# ## cannot run clients on it. Since NeAT relies on WSDL clients, it is
# ## impossible to have neat running on and Ubuntu server. The
# ## installation is however possible, since the stubb can be generated
# ## on rsat-tagc.univ-mrs.fr.  I have no idea how we did to install
# ## SOAP::WSDL there. In any case, the 
# ##
# ## Solution proposed here: http://stackoverflow.com/questions/3489642/dependency-problem-of-perl-cpan-modules
# ## Not sure it works by its own, but cannot harm.
# cpan
# ## At the cpan prompt, type the following
# install Module::Build
# install Module::Build::Compat
# install CPAN ## This takesa HUGE time. I answer all questions by the default answer (simply type the Enter key)
# upgrade ## Takes a HUGE time, since all packages are apparently re-tested
# quit
##
## DONE: tricky stuff
################################################################

################################################################
################       RSAT installation        ################        TO BE REMOVED ? TO BE UPDATED
################################################################

## Create a specific user for RSAT. The user is named rsat
sudo adduser rsat
## Full Name: Regulatory Sequence Analysis Tools admin

## Grant sudoer privileges to the rsat user (will be more convenient for
## installing Perl modules, software tools, etc)
sudo visudo
## then add the following line below "User privilege specification"
# rsat    ALL=(ALL:ALL) ALL

################################################################
## Download RSAT distribution

## Note: the git distribution requires an account at the ENS git
## server, which is currently only possible for RSAT developing team.
## In the near future, we may use git also for the end-user
## distribution.

## RSAT installation is done under the rsat login.
##
## We retrieve it in the home folder of the RSAT user, because RSAT ##### mot de passe git ? clé publique ??
## user has no write authorization on /bio. In a second time we do a
## sudo mv to place the rsat folder in /bio.
# su - rsat
# cd ${HOME}
# git clone git@depot.biologie.ens.fr:rsat ## euh password ?!
# git config --global user.mail rsat@rsat-vm-2015-02
# git config --global user.name "RSAT admin on RSAT-VM"
# 

## Define an environment variable with the RSAT_HOME directory
## (will be used later to configure RSAT)
export INSTALL_ROOT_DIR=/bio
export RSAT_HOME=/bio/rsat

## Move the rsat distribution to the RSAT_HOME directory
# sudo mv ${HOME}/rsat /bio/rsat
# ln -fs /bio/rsat ${HOME}/rsat

## For users who don't have an account on the RSAT git server, the
## code can be downloaded as a tar archive from the Web site.
##
## This is however less convenient than using the git clone, which
## greatly facilitates updates.
#
# ! voir plus haut
# export RSAT_DISTRIB=rsat_2014-08-22.tar.gz
# export RSAT_DISTRIB_URL=http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/${RSAT_DISTRIB}
#
export RSAT_DISTRIB=rsat_2015-06-06.tar.gz
export RSAT_DISTRIB_URL=http://pedagogix-tagc.univ-mrs.fr/download_rsat/${RSAT_DISTRIB}

cd /bio
sudo mkdir -p /bio/rsat
sudo wget ${RSAT_DISTRIB_URL}
sudo tar -xpzf ${RSAT_DISTRIB}
sudo rm -f   ${RSAT_DISTRIB} ## To free space
sudo df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_downloaded.txt
cd ~; ln -fs /bio/rsat rsat ## -> link to self ? error: ln: failed to access ���rsat���: Too many levels of symbolic links
                        

## Metabolic pathway tools installation
##
## TO BE DONE LATER
## wget http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/metabolic-tools_20110408.tar.gz

# ## DONE : obtained RSAT distribution
# ################################################################

## Run the configuration script, to specify the environment variables.
cd /bio/rsat
sudo perl perl-scripts/configure_rsat.pl ## what params ? 

## Parameters to change
##   rsat_site            rg-vm
##   rsat_server_admin    my.email@inserm.fr

## Load the (updated) RSAT environment variables
source RSAT_config.bashrc

## Check that the RSAT environment variable has been properly configured
echo ${RSAT}

######################################################################################
## Grant sudoer privileges to the rsat user (will be more convenient for
## installing Perl modules, software tools, etc)
sudo visudo
## then add the following line below "User privilege specification"
# rg    ALL=(ALL:ALL) ALL
###################################################################################### TO BE UPDATED



## Initialise RSAT folders
make -f makefiles/init_rsat.mk init

################################################################
## For the next operations, we need to be su

## Since we became rsat user with "sudo", a simple exit brings us back
## to the root user
exit

## Load RSAT bashrc file
cd /bio/rsat
source RSAT_config.bashrc

## Check that the root has well loaded the RSAT configuration
echo $RSAT

################################################################
## Previous way to specify bashrc parameters, via
## /etc/bash_completion.d/. I change it (2014-09-23) because it does
## not allow to run remote commands via ssh (/etc/bash_completion.d is
## apparently only loaded in interactive mode).
## 
## Link the RSAT bash configuration file to a directory where files
## are loaded by each user at each login. Each user will then
## automatically load the RSAT configuration file when opening a bash
## session.


rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/                      #### ??
## ln -fs /bio/rsat/RSAT_config.bashrc /etc/bash_completion.d/

#emacs -nw /etc/bash.bashrc




################################################################
## Installation of Perl modules required for RSAT
################################################################
## Notes
##
## 1) limxml2-dev is required to compile the Perl module XML::LibXML
##        sudo apt-get install limxml2-dev 
##
## 2) For some modules, installation failed until I used "force"
##	 force install SOAP::WSDL
##	 force install SOAP::Lite
##
## 3) Problem of dependency when installint REST::Client, even with "force".
##	 MCRAWFOR/REST-Client-271.tar.gz              : make_test FAILED but failure ignored because 'force' in effect
##	 NANIS/Crypt-SSLeay-0.64.tar.gz               : make NO
## To solve it I found this
## 	sudo apt-get install libnet-ssleay-perl
##	sudo apt-get install libcrypt-ssleay-perl
##	sudo apt-get install libssl-dev


## The first time we use cpan, we apparently need to open it manually
## in order to configure and update it.
cpan                              ### params -> default
install YAML
## I am not sure, but I think that this command is useful to properly install the subsequent packages.
install CPAN 
reload cpan
## then type "quit"

## Set working directory to RSAT
cd $RSAT

## Get the list of Perl modules to be installed
make -f makefiles/install_rsat.mk  perl_modules_list

## Check which Perl modules are already installed
make -f makefiles/install_rsat.mk perl_modules_check
## The result file name will be displayed at the end of the tests

## Check the result of perl modules test
more check_perl_modules_eval.txt
## On Ubuntu 14.04, Object::InsideOut has status "Fail" but there is
## apparently no problem

## Check missing modules
grep Fail  check_perl_modules_eval.txt

## Identify Perl modules that were not OK after the ubuntu package installation
grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;'
MISSING_PERL_MODULES=`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | cut -f 2 | xargs`
echo "Missing Perl modules:     ${MISSING_PERL_MODULES}"


## Beware: the _noprompt suffix is
## optional. It has the advantage to avoid for the admin to confirm
## each installation step, but the counterpart is that errors may be
## overlooked.
make -f makefiles/install_rsat.mk perl_modules_install PERL_MODULES="${MISSING_PERL_MODULES}"

## Check if all required Perl modules have now been correctly installed
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt
## Note: Object::InsideOut should be ignored for this test, because it
## always display "Fail", whereas it is OK during installation.

## Note: I had to force installation for the some modules, because
## there seem to be some circular dependencies.
grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v "Object::InsideOut"
make -f makefiles/install_rsat.mk perl_modules_install_by_force PERL_MODULES_TO_FORCE="`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v Object::InsideOut| cut -f 2 | xargs`"


## Last check for Perl modules. 
## If some of them still fail (except Object::InsideOut), manual intervention will be required.
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt

## Measure remaining disk space
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_perl_modules_installed.txt
grep ${DEVICE} /bio/install_logs/df_*.txt


################################################################
## Configure RSAT web server

## Edit the file to replace [RSAT_PARENT_FOLDER] byt the parent directory
## of the rsat directory.
cd ${RSAT}
sudo rsync -ruptvl RSAT_config.conf /etc/apache2/sites-enabled/rsat.conf

## OPTIONAL: since I am using this to install a virtual machine whose
## only function will be to host the RSAT server, I replace the normal
## default web folder by RSAT web folder. 
##
emacs -nw /etc/apache2/sites-available/000-default.conf
## Comment the line with the default document root (should appear as
## such in the original config):
##        DocumentRoot /var/www/html                                                                            
## And write the following line:
##        DocumentRoot /bio/rsat/public_html
## The server will now immediately display RSAT home page when you
## type its IP address.

apache2ctl restart

## You should now test the access to the RSAT Web server, whose URL is
## in the environment variable RSAT_WWW
echo $RSAT_WWW

## If the value is "auto", get the URL as follows
export IP=`ifconfig eth0 | awk '/inet /{print $2}' | cut -f2 -d':'`
# export IP=192.168.56.101
echo ${IP}
export RSAT_WWW=http://${IP}/rsat/
echo $RSAT_WWW

################################################################
## Next steps require to be done as rsat administrator user

## compile RSAT programs written in C
su - rsat
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all
export INSTALL_ROOT_DIR=/bio/
sudo df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_app_compiled.txt

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!  BUG    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!! I HAVE A PROBLEM TO COMPILE KWALKS. SHOULD BE CHECKED !!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Install some third-party programs required by some RSAT scripts.
make -f makefiles/install_software.mk install_ext_apps
sudo df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_extapp_installed.txt

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!      ONLY FOR THE IFB CLOUD    !!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##
## replace the data directory by a link to
## a separate disk containing all RSAT data.
export RSAT_DATA_DIR=/root/mydisk/rsat_data
cd ${RSAT}/public_html
mv data/* ${RSAT_DATA_DIR}/
mv data/.htaccess ${RSAT_DATA_DIR}/
rmdir data
ln -s ${RSAT_DATA_DIR} data
cd $RSAT

## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae \
 -org Escherichia_coli_K_12_substr__MG1655_uid57779

## Optionally, install some pluricellular model organisms
# download-organism -v 1 -org Drosophila_melanogaster
# download-organism -v 1 -org Caenorhabditis_elegans
# download-organism -v 1 -org Arabidopsis_thaliana

## Get the list of organisms supported on your computer.
supported-organisms

sudo df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_organism_installed.txt

################################################################
## Install selected R librairies, required for some RSAT scripts
################################################################

whoami
## Should return "rsat"

## Installation of R packages
cd $RSAT; make -f makefiles/install_rsat.mk install_r_packages

## More convenient: the following command does the same (install R
## packages) + compile the C programs
cd $RSAT; make -f makefiles/install_rsat.mk update

# ## The command R CMD INSTALL apparently does not work at this stage.
# ##	root@rsat-tagc:/workspace/rsat# R CMD INSTALL reshape
# ##	Warning: invalid package 'reshape'
# ##	Error: ERROR: no packages specified
# cd $RSAT; make -f makefiles/install_rsat.mk  r_modules_list 
# ### I install them from the R interface. This should be revised to
# ### make it from the bash, but I need to see how to specify the CRAN
# ### server from the command line (for the time being, I run R and the
# ### programm asks me to specify my preferred CRAN repository the first
# ### time I install packages).
# R
# ## At the R prompt, type the following R commands.
# ##
# ## Beware, the first installation of bioconductor may take a while,
# ## because there are many packages to install
# ##
# ## Note: since this is the first time you install R packages on this
# ## VM, you need to choose a RCRAN server nearby to your site.
# install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras", "dendextend"))
# source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
# quit()
# ## At prompt "Save workspace image? [y/n/c]:", answer "n"

## Check remaining disk space
df -m > /bio/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_R_packages_installed.txt

################################################################
## At this stage you can already check some simple RSAT command 

## Test a simple Perl script that does not require for organisms to be
## installed.
which random-seq
random-seq -l 100

## Test a simple python script that does not require organisms to be
## installed.
random-motif -l 10 -c 0.90

################
## Test some external programs

## vmatch (used in purge-sequence)
random-seq -l 100 | purge-sequence

## Check that seqlogo is installed
which seqlogo
seqlogo

## Check that weblogo 3 is installed
which weblogo
weblogo --help

## ghostscript
which gs
gs --version

## Check tat the model genomes have been correctly installed
## A simple and quick test: retrieve all the start codons and count
## oligonucleotide frequencies (most should be ATG).
retrieve-seq -org Saccharomyces_cerevisiae -all -from 0 -to +2 \
    | oligo-analysis -l 3 -1str -return occ,freq -sort

################################################################
## Configure the SOAP/WSDL Web services

## Check the URL of the web services (RSAT_WS). By default, the server
## addresses the WS requests to itself (http://localhost/rsat) because
## web services are used for multi-tierd architecture of some Web
## tools (retrieve-ensembl-seq, NeAT).
cd $RSAT
#echo $RSAT_WS

## Get tehe current IP address
export IP=`/sbin/ifconfig eth0 | awk '/inet /{print $2}' | cut -f2 -d':'`
# export IP=192.168.56.101
echo ${IP}
export  RSAT_WS=http://${IP}/rsat/
## Initialize the Web services stub. 
make -f makefiles/init_rsat.mk ws_init RSAT_WS=${RSAT_WS}

## After this, re-generate the web services stubb, with the following
## command.
make -f makefiles/init_rsat.mk ws_stub RSAT_WS=${RSAT_WS}

## Test the local web services
make -f makefiles/init_rsat.mk ws_stub_test

## Test RSAT Web services (local and remote) without using the
## SOAP/WSDL stubb (direct parsing of the remote WSDL file)
make -f makefiles/init_rsat.mk ws_nostub_test

## Test the program supported-organisms-server, which relies on Web
## services without stub
supported-organisms-server -url ${RSAT_WS} | wc
supported-organisms-server -url http://localhost/rsat/ | wc
supported-organisms-server -url http://rsat-tagc.univ-mrs.fr/ | wc


################################################################
## tests on the Web site

## Run the demo of the following tools
##
## - retrieve-seq to check the access to local genomes (at least
##   Saccharomyces cerevisiae)
##
## - feature-map to check the GD library
##
## - retrieve-ensembl-seq to check the interface to Ensembl
##
## - fetch-sequences to ceck the interface to UCSC
##
## - some NeAT tools (they rely on web services)
##
## - peak-motifs because it mobilises half of the RSAT tools -> a good
##   control for the overall installation.
##
## - footprint-discovery to check the tools depending on homology
##   tables (blast tables).


################################################################
## Install the cluster management system (torque, qsub, ...)

## Check the number of core (processors)
grep ^processor /proc/cpuinfo

## Check RAM
grep MemTotal /proc/meminfo


################################################################
################ Install Sun Grid Engine (SGE) job scheduler
################################################################


## Beware, before installing the grid engine we need to modify
## manually tjhe file /etc/hosts
emacs -nw /etc/hosts
## Initial config (problematic) 
##    127.0.0.1       localhost       rsat-vm-2015-02
##    127.0.1.1      rsat-vm-2015-02
## Config to obtain: 
##    127.0.0.1       localhost       rsat-vm-2015-02
##    #127.0.1.1      rsat-vm-2015-02
apt-get install --quiet --assume-yes gridengine-client
apt-get install --quiet --assume-yes gridengine-exec
apt-get install --quiet --assume-yes gridengine-master
apt-get install --quiet --assume-yes gridengine-qmon 

qconf -aq default  ## aggregate a new queue called "default"
qconf -mq default  ## modify the queue "default"
qconf -as localhost ## aggregate the localhost tho the list of submitters
## -> set the following values
## hostlist              localhost

## Take all default parameters BUT For the SGE master parameter, type
## localhost (it must be the hostname)

## Test that jobs canbe sent to the job scheduler


################################################################
########################     OPTIONAL     ######################
################################################################


## Install some software tools for NGS analysis
cd ${RSAT}
make -f makefiles/install_software.mk install_meme

################################################################
## Ganglia: tool to monitor a cluster (or single machine)
## https://www.digitalocean.com/community/tutorials/introduction-to-ganglia-on-ubuntu-14-04
sudo apt-get install -y ganglia-monitor rrdtool gmetad ganglia-webfrontend
sudo cp /etc/ganglia-webfrontend/apache.conf /etc/apache2/sites-enabled/ganglia.conf
sudo apachectl restart


