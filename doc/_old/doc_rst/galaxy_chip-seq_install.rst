Downloading Galaxy code
=======================

We followed the instructions from the Galaxy Web site:

-  https://wiki.galaxyproject.org/Admin/GetGalaxy

\`\`\`{r eval=FALSE} ## get a git clone of galaxy git clone
https://github.com/galaxyproject/galaxy/ cd galaxy ## Go th the galaxy
directory

Check out the master branch, recommended for production server
==============================================================

| git checkout -b master origin/master
| git pull ## Just in case, we are already up-to-date \`\`\`

Configure the Galaxy server (and get python modules if required)
================================================================

We first edit the config file to chooe a specific port for Galaxy

``{r eval=FALSE} cp config/galaxy.ini.sample config/galaxy.ini``

We then edit this file by setting the port to 8082, because our 8080 is
already used for other purposes.

We performed the following modifications.

admin\_users=admin1@address.fr,admin2@univbazar.fr,admin3@gmail.com port
= 8082 # The port on which to listen. host = 0.0.0.0 ## To enable access
over the network allow\_user\_deletion = True

Configuring the Apache server on RSAT
=====================================

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
==========================

``{r eval=FALSE} sh run.sh``

On our internal network, the server becomes available at the address:

http://192.168.1.6:8082

Registrating
============

-  open a connection to the Galaxy server
-  In the Galaxy menu, run the command **User -> Register**. Enter the
   same email address as you declared as admin users.

Install Galaxy modules
======================
