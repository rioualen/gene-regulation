
# Creating an NGS VM under Virtualbox


## To do

1. Create a separate virtual disk for the data, and explain how to mount it on the virtual machine.

## Requirements

1. Virtualbox software can be downloaded from
    [https://www.virtualbox.org/](https://www.virtualbox.org/)

2. The VM is built under Ubuntu 14.04 (latest long-term supported version)

    [http://releases.ubuntu.com/14.04/ubuntu-14.04.2-desktop-amd64.iso](http://releases.ubuntu.com/14.04/ubuntu-14.04.2-desktop-amd64.iso)

------------------------------------------


## Creation of the virtual machine

1. Open VirtualBox

2. Click on the **New** button. 

3. Parameters

- Name and operating system

    - Name: reg-genomics-vm-2015-07
    - Type: Linux
    - Version: Ubuntu (64 bits)
    
- Memory size:	2048 Mb (this can be modified afterwards). 

- Hard drive:	*Create a virtual hard drive now*. 

- Hard drive file type: *VMDK* (Virtual Machine Disk). 

    I chose this option because it ensures a wider compatibility with other OS and Virtual Machine management systems. 

	   Another potential advantage of VMDK is that it enables to split virtualdisks in files <=2Gb, which is convenient to store them on FAT partitions. 

- Storage on physical hard drive
    - Select *Dynamically allocated*
    - Activate the option *Split into files less than 2Gb*, which allows to store the VM on FAT partitions for Windows host machines.


- File location and size

    - name of virtual hard drive file:    reg-genomics-hdd
    - max size of virtual hard drive:	    18GB

        Beware: the size has to be larger than you real needs, because Ubuntu will automatically create some big  partitions (/dev and /run/shm). Since we allow dynamical memory allocation, it is fine to set the max size to 18Gb, which will be grossly assigned as follows:

        - max 14Gb for the root partition /
        - 2Gb for /dev
	      - 2Gb for /run/user




At this stage, the VM has been created and needs to be configured before installing the operating system. 

---------------------------------------


## General configuration of the network for your VirtualBox program

Before configuring the virtual machine, we need to tell VirtualBox how it will enable your local virtual machines to interact with their host (the operating system of the machine on which the VM is running).

1. Open **VirtualBox Preferences**

2. Open the tab *Network* > *Host-only Networks*
     - click on the "+" icon
     - this creates a network vboxnet0. Select this network, click on the screw driver icon (*edit host-only network*), and set the following options:

      - *Adapter* tab
          - IPv4 Address: 192.168.56.1
          - IPv4 Network Mask: 255.255.255.0
          - IPv6 Adress: blank
          - IPv6 Network Mask Length: 0

      - *DHCP Server* tab
          - Check *Enable Server*
          - *Server Address:* 192.168.56.100
          - *Server Mask:* 255.255.255.0
          - *Lower Address Bound:* 192.168.56.101
          - *Upper Address Bound:* 192.168.56.254

## Configuration of the virtual machine

In the VirtualBox main window, select the newly created virtual machine, and click on the **Settings** button. We will successively select the different tabs and configure them. 

### Network

VirtualBox offers many alternative ways to configure network communications between the virtual machine, the host machine, and the external network. 

To get more information about network settings: 
    - VirtualBox manual page: [https://www.virtualbox.org/manual/ch06.html](https://www.virtualbox.org/manual/ch06.html)).
    - An excellent tutorial: [http://christophermaier.name/blog/2010/09/01/host-only-networking-with-virtualbox](http://christophermaier.name/blog/2010/09/01/host-only-networking-with-virtualbox)

We present here one possible way to configure your Virtual machine, but this should be adapted to the particular security/flexibility requirements of the network where the maching has to run. 


In the VM settings, select tne *Network* tab. VirtualBox enables you to specify several adapters, each corresponding to one separate network acces (e.g. using an ethernet card + wi-fi connection).

- click on the tab *Adapter 1*,
- check *Enable Network Adapter*
- Attached to: *Host-only Adapter*
- Name: *vboxnet0* (this network must have been created beforehand)


#### Additional networks (optional)


**BEWARE!** If you activate his adaptor your VM will be visible from external computers in the same network as the host  machine.

See information at: [https://www.virtualbox.org/manual/ch06.html#network_bridged](https://www.virtualbox.org/manual/ch06.html#network_bridged)

- click on the tab *Adapter 2*,
    - check *Enable Network Adapter*
    - Attached to : *Bridged Adapter*
    - Name: choose an option corresponding to the actual internet connection of the host machine (e.g. ethernet cable, Wi-Fi, ...).

Note: you can also activate a third adapter, for example to support both an ethernet card (Adapter 2) and a Wi-Fi connection (Adapter 3). This is particularly interesting for laptop computers. 



================================================================
Ubuntu operating system installation
================================================================

Desktop version
---------------

The desktop version of Ubuntu presents the advantage of letting users
work directly in the graphical environment of the Ubuntu virtual
machine, which is convenient to open non-text result files for
visualization (web pages, images) and to analyze the results in other
graphical packages (R, openOffice, CytoScape, ...).

Ubuntu install disk downloaded from
       http://releases.ubuntu.com/14.04.1/ubuntu-14.04.1-desktop-amd64.iso

Started the default installation. Steps:
	- click the button: install Ubuntu

Preparing to install Ubuntu
	  - leave all default parameters and click GO

Installation type
	- (leave the default) Erase disk and install Ubuntu, click "Install Ubuntu"

Where are you (automatic)
      - 2015-02 I was in Mexico

Keyboard layout

      - I have a French - French Macintosh, so I use it for
        installation, but I will then add an English keybord.

NOTE CLAVIER FRANÇAIS: for the French Macintosh keyboard, I spent a
lot of time to find the solution to use the 3-component keys
(e.g. alt-shit-L for the pipe character |). The simplest solution: use
the alt key on the RIGHT side of the space bar.

Who are you ?
    Your name			RSAT admin
    Your computer's name	rsat-vm-2015-02
    Pick a username 		rsat

    Change since 2015-02: 
    	   I activate the option Log in automatically
	   This will greatly facilitate the use of the Destkop version

# I should think about the possibility to activate "Log in
# automatically", in order to facilitate the distribution without
# providing a password to the users.



================================================================
Screen resolution (for the Desktop version)
================================================================

At first log in, screen resolution is restricted to 640x480 pixels. 

I found the solution here:

http://www.juanrubio.me/2014/02/ubuntu-trusty-virtualbox-guest-additions/

I run the following commands: 

  sudo bash
  apt-get install virtualbox-guest-additions-iso
  software-properties-gtk --open-tab=4

In the dialog box, 
   - check the option "Using x86 Virtualization ..."
   - click "Apply changes" 
   - click "Close"
Reboot


================================================================
====
==== Server version
====
================================================================

The server version is less comfortable than the desktop version, but
it should be more economic in terms of hard drive (and CPU/RAM for
graphical resources ?).

Language  English

Home page: choose "Install Ubuntu Server"

Language (again): English

Country, territory or area: other > Europe > France

Country to base the default locale setting: United Kingdom 
	(there is apparently no setting for France !)

Detect keyboard layout: No > French > French - French (Macintosh)

Your system has multiple network interfaces (I guess this corresponds
tothe 3 adapters specified in VirtualBox preferences).
      Primary network interface: eth0

host name: rsat-vb-ub14s

Full name of new user:	RSAT admin

Username for your account: rsat

Password for new user: the classical (only for RSAT team members)

Encrypt home directory ? No

Time: Europe/Paris

Partitioning method: (default) Guided - use entire disk and set up LVM

     Next partitioning options: I say "continue" to all with default values

 ... Installing the sysem ... (takes some time)

Configure the package manager
	  Leave the HTTP proxy blank
  -> ... configuring apt (takes some time)

Configuration taskse1
   Install security updates automatically

Software selection : I just install the OpenSSH server, in order to be
able connecting, but I don't install anything else since I want to use
the RSAT bash file with apt-get packages etc.
    -> ... Select and install software (takes some time)

Install the GRUB boot loader
	Yes

================================================================
Login in in the VM
================================================================

At this stage, the machine is able to boot and get a dynamic IP from
the network. I login as rsat user and type
    ifconfig
to get the IP address. I then open a terminal on the host machine, and run
   ssh rsat@[IP_OF_MY_VM]
so I can run the rest of the installation in my familiar environment
(terminal, keyboard, ...)


ssh rsat[myvm]

## Transfer my screenrc preferences
rsync -ruptvl rsat@rsat.ulb.ac.be:.screenrc .

## Transfer RSAT ssh parameters + open an agent that will manage my
## password for future ssh communications.
rsync -ruptvl rsat@rsat.ulb.ac.be:.ssh .
ssh-agent > ~/agent
source ~/agent
ssh-add

================================================================
Network specification
================================================================


## We need to install the ssh package in order to enable connection
## from the host terminal
apt-get install --quiet --assume-yes ssh


## We can now access the VM by ssh to its static IP address
ssh rsat@192.168.56.101

################################################################
### ### NOT NECESSARY ANYMORE

### ## I need emacs to edit the network configuration file
### sudo bash
### df -h
### apt-get install  --quiet --assume-yes  emacs

### ## Since this is the first installation of an apt-get package, it
### ## takes some time to install dependencies.
### df -h

### emacs -nw /etc/network/interfaces

### ## I edit the file /etc/network/interfaces

### ## Host-only network interface for the VirtualBox VM
### auto eth0
### iface eth0 inet static
### address 192.168.56.114
### netmask 255.255.255.0
### network 192.168.56.0
### broadcast 192.168.56.255
### gateway 192.168.56.1

### ## Let the second interface automatically ask an IP address to the
### ## DHCP server (this was by default the first interface eth0, but we
### ## want to impose the static address as first interface).
### autho eth1
### iface eth1 inet dhcp

### ## Re-activate the interface, this seems to be the best way withb ubuntu 14.04
### ifdown eth1 ## Stop the eth1 interface
### ifup eth1 ## Start the eth1 interface with its new configuration
### ifconfig ## Check that eth1 is specified as expected

### ## This should display something like this
### ## ...
### ## eth1      Link encap:Ethernet  HWaddr 08:00:27:9e:3d:cd  
### ##           inet addr:192.168.56.114  Bcast:192.168.56.255  Mask:255.255.255.0
### ## ...

### ## reboot
### ## (Note: I am not sure if we still need to reboot here, to be checked)

================================================================
Install ubuntu packages
================================================================

## cf script
##   $RSAT/doc/howto/install_rsat_ubuntu14.04.bash


================================================================
Useful commands to manage VirtualBox
================================================================

List running VMs
  VBoxManage list runningvms

Get the IDs of the virtual machines
  VBoxManage list runningvms | awk -F"[{}]" '{print $2}'

For convenience, we store the IP of the machine of interest in an
environment variable VMID.
  VMID=`VBoxManage list runningvms | awk -F"[{}]" '{print $2}'`

Get properties of the VM
  VBoxManage guestproperty enumerate ${VMID}

List information about a given virtual machine
  VBoxManage showvminfo ${VMID}

Note: that command fives the MAC address, but no IP address.
  VBoxManage showvminfo ${VMID} | grep '^NIC'

================================================================
Convenient settings (not necessary)
================================================================


1) Clipboard sharing between guest and host
-------------------------------------------
For the desktop version of Ubuntu, it is convenient to enable
copy-paste between the guest and the host.

- select the VM
- In the VirtualBox menu, select "Machine -> Settings"
- In panel "General", select the tab "Advanced".
     set "Shared clipboard" to "Bidirectional"

================================================================
Create a generic user for the virtual machine
================================================================


################################################################
## Create a user for the virtual machine
##
## This VM user is separate from the rsat user, which only serves to
## manage the RSAT software suite and related packages.
##
## For the sake of security, we force this user to change password at
## first login

## First delete this user (in case it was previously defined)
##  sudo userdel --remove vmuser

sudo useradd --password `openssl passwd -1 -salt xyz tochng`\
    --home /home/vmuser \
    --create-home \
    --shell /bin/bash \
    --comment "VM user" \
    vmuser

## Force vmuser to change password at first login
sudo chage -d 0 vmuser

## Add sudoer rights to vmuser
sudo chmod 644 /etc/sudoers
sudo emacs -nw /etc/sudoers
## Find the following line
##     # User privilege specification
##     root    ALL=(ALL:ALL) ALL
## Below it, add the following line:
##     rsat  ALL=(ALL:ALL) ALL
##     vmuser  ALL=(ALL:ALL) ALL


================================================================
PROBLEMS TO BE FIXED
================================================================

1) Ethernet: optimize solution for a simple setting in the classroom.

2) Problem of disk occupancy: Ubuntu shared memory occupies 2Gb. This
   can be modified as explained here.
	http://www.cyberciti.biz/tips/what-is-devshm-and-its-practical-usage.html

4) For the desktop version, I should try to use the "Advanced
   installation" in order to use less disk without loosing too much
   confort. For example, I could inactivate the support for all the
   languages that are installed by default (downloading language
   packs).


################################################################
##
## Customizing your instance of the VM
##
################################################################


## Set up time zone, date and time (source:
## https://help.ubuntu.com/community/UbuntuTime).
sudo dpkg-reconfigure tzdata


################################################################
## To set up the keyboard (may vary between users who download the
## VirtualBox VM)
sudo dpkg-reconfigure console-data
## (for my laptop, mac / Unknown / French /Standard / New)


################################################################
## For the Virtualbox VM: mount the virtual disk rsat_data



================================================================
Export appliance
================================================================

Once the Virtual Machine is working fine, it can be exported to an
appliance.

In VirtualBox, open "File -> Export Appliance ...", and select the VM.

   Save as: rsat-vm-2015-02_[YYYY-MM-DD].ova
   Files of Type: Open Virtualization Format Archive (.ova)

Storage Settings
	Format: OVF 1.0
	v Write Manifest File


Appliance Settings
	Name		rsat-vb-ub14s ***OR***  rsat-vm-2015-02
	Product		Regulatory Sequences Analysis Tools
	Product-URL	http://rsat.eu/
	Vendor		Jacques van Helden
	Vendor URL	http://jacques.van-helden.perso.luminy.univ-amu.fr/
	Version		2014-10-10
	Description	Regulatory Sequence Analysis Tools (RSAT, http://rsat.eu/), installed on an Ubuntu 14.04 Virtual Machine. 
	License		Free of use for academic users, non-commercial and non-military usage. 



################################################################
## Prepare the VM + accompanying material on a USB key, to ensure
## distribution in the teaching room.
##
## We use ScanDisk 64Gb USB3.0 keys, which have a high transfer rate
## (190Mbps write / 250Mbps read)

SOURCE_DIR=/no_backup/VirtualBox_VMs/appliances
VERSION=2014-09-05a

## Specification of the target location
KEY_NB=02
TARGET_DIR=/Volumes/RSAT-VM_${KEY_NB}/
time rsync -ruptvl ${SOURCE_DIR} ${TARGET_DIR}

################
# OPTIONAL: split the VM into smaller files. This can be useful for two purposes.
#
# 1) To transfer files to a FAT-formatted drive for distribution (FAT
#    cannot hold files >4Gb). Note that the final installation drive
#    will have to support files >4Gb for the server version of the VM.
# 2) To facilitate synchronization between local and remote machines.

SPLIT_SIZE=3999m
mkdir -p ${TARGET_DIR}/appliances
time split -b ${SPLIT_SIZE} ${SOURCE_DIR}/rsat-vm-2015-02_${VERSION}.ova ${TARGET_DIR}/appliances/rsat-vm-2015-02_${VERSION}.ova.split${SPLIT_SIZE}_
################


