# pcaExplorer - Interactive exploration of Principal Components of Samples and Genes in RNA-seq data - *Server Edition*

`pcaExplorer` is a Bioconductor package containing a Shiny application for analyzing expression data
in different conditions and experimental factors.

It is a general-purpose interactive companion tool for RNA-seq analysis, which guides the user in 
exploring the Principal Components of the data under inspection. `pcaExplorer` provides tools and 
functionality to detect outlier samples, genes that show particular patterns, and additionally 
provides a functional interpretation of the principal components for further quality assessment and 
hypothesis generation on the input data.

This repository contains the information and the files for setting up a Shiny Server instance
where `pcaExplorer` runs for you.

The code contained here is also released under the MIT License.

More information on `pcaExplorer` can be found at https://github.com/federicomarini/pcaExplorer - as 
well at http://bioconductor.org/packages/pcaExplorer/

## Setting up the Shiny Server

Mainly refer to this excellent reference: http://docs.rstudio.com/shiny-server/

In my case, I tested the deployment on a Ubuntu-based server. If this is your case too, please go to
http://docs.rstudio.com/shiny-server/#ubuntu-12.04

If you are already familiar with Shiny Server, you can consider going directly here: 
http://docs.rstudio.com/shiny-server/#quick-start

Follow the indications in the main reference guide to see whether you have all packages installed, 
and visit the address which will refer to the address of your server 
(http://your.server.address:3838/sample-apps/hello/). 


## Setting up `pcaExplorer`


You need to:

- install `pcaExplorer`
- setup `pcaExplorer` on the server

### Installing `pcaExplorer`

You need first to install Bioconductor

```
sudo su - -c "R -e \"source('http://bioconductor.org/biocLite.R')\""
```

Once that is done, you can install `pcaExplorer` from GitHub with

```
# install devtools first
sudo su - -c "R -e \"BiocInstaller::biocLite('devtools')\""
# and then pcaExplorer
sudo su - -c "R -e \"BiocInstaller::biocLite('federicomarini/pcaExplorer')\""
```

`pcaExplorer` is also on Bioconductor, you can install with

```
sudo su - -c "R -e \"BiocInstaller::biocLite('pcaExplorer')\""
```


### Setup `pcaExplorer` on the server

You essentially need to do two things:

- in `/srv/shiny-server`, 

```
git clone https://github.com/federicomarini/pcaExplorer_serveredition.git
```

This will clone the repository with all the files required for running

- you need to edit `/etc/shiny-server/shiny-server.conf`. To do this, use any text editor (e.g. `vi`),
and add the following lines 

```
server {
  listen 3838;

  #...

  ### FROM HERE
  location /pcaExplorer {
    # Run this location in 'app_dir' mode, which will host a single Shiny
    # Application available at '/srv/shiny-server/myApp'
    app_dir /srv/shiny-server/pcaExplorer;

    # Log all Shiny output to files in this directory
    log_dir /var/log/shiny-server/pcaExplorer;

    # When a user visits the base URL rather than a particular application,
    # an index of the applications available in this directory will be shown.
    directory_index off;

    # recommended, to
    app_init_timeout 250;
  }
  ### TO HERE

  # ...
}
```

This is basically as suggested in the Shiny Server reference, with the only exception of the 
`app_init_timeout` value: set this a little higher, so that you avoid that the service stops too often.
Credits go to Mark Robinson for suggesting this parameter setup!


You might need to restart the Shiny Server (`systemctl restart shiny-server` on Ubuntu). Then you should 
be good to go!

Visit http://your.server.address:3838/pcaExplorer for your personal running instance.


## `pcaExplorer` at the IMBEI

http://shiny.imbei.uni-mainz.de:3838/pcaExplorer is the address for the public instance of `pcaExplorer` on the 
Shiny Server I manage at the Institute of Medical Biostatistics, Epidemiology and Informatics.

The relevant setup of that machine is the following (in case you are asking what specs you need):

- 4 cores
- 8 GB RAM

```
lsb_release -a

	No LSB modules are available.
	Distributor ID:	Ubuntu
	Description:	Ubuntu 16.04.2 LTS
	Release:	16.04
	Codename:	xenial
```






## Contact

For additional details regarding the functions of `pcaExplorer`, please consult the package documentation or 
write an email to marinif@uni-mainz.de. 

For questions regarding the deployment on your Shiny Server, please use 
https://github.com/federicomarini/pcaExplorer_serveredition/issues


### Bug reports/Issues/New features 

Please use https://github.com/federicomarini/pcaExplorer/issues for reporting bugs, issues or for 
suggesting new features to be implemented.
