---
output: html_document
---

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


### Storage

Click on the **Empty** disc icon in the storage tree.  Select the disc icon on the right and fetch the downloaded \*.iso image  (see **Requirements**)
Click on *OK*

You can now start the VM. 

## Operating system installation

* Welcome 
    - check the language settings and click on *Install Ubuntu*.

* Preparing to install Ubuntu
	  - leave all default parameters and click *Continue*.

* Installation type
	  - (leave the default) Erase disk and install Ubuntu, click *Install Now*.

* Where are you (automatic)
      - Paris

* Keyboard layout

      - French - French

* Who are you ?
      - Your name			        Regulatory Genomics
      - Your computer's name	reg-genomics
      - Pick a username 		  rg
      - Choose a password     reg-genomics
      - (Activate the option Log in automatically)
      
Restart once installation is completed. 

Once on the desktop, go to the VM menu: select *Devices* then *Install Guest Additions CD image*.
Run it. 

The VirtualBox Guest Additions will provide closer integration between host and guest and improve the interactive performance of guest systems.
Reboot again to see the new display.


## Connection to VM through the host machine. 

In the VM terminal: 

```
apt-get install --quiet --assume-yes ssh ## install ssh
ifconfig ## get the IP
```

In the host's terminal: 

```
ssh rsat@192.168.56.101 ## access the VM by ssh using its IP address
```

```
rsync -ruptvl {me}@{myIP}:.ssh .
ssh-agent > ~/agent
source ~/agent
ssh-add
```
#### OK til here


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



#### Customizing your instance of the VM

Set up time zone, date and time (source:
https://help.ubuntu.com/community/UbuntuTime).
sudo dpkg-reconfigure tzdata


To set up the keyboard (may vary between users who download the VirtualBox VM)
sudo dpkg-reconfigure console-data
(for my laptop, mac / Unknown / French /Standard / New)

For the Virtualbox VM: mount the virtual disk rsat_data



### Export appliance

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


