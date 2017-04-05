Gene-regulation library
================================================================

The library contains a variety of files, including scripts and 
configuration files. 

You can find a description of these hereafter.

Snakemake files (snakefiles)
----------------------------------------------------------------

Snakefiles are based on the scripting language Python 3, and use a specific syntax.

For organization purpose, we have distinguished two types 
of Snakefiles: 

* Rules are typically "bricks" to build workflows with. 
Eache rule corresponds to a specific operation.
however one wants.

* Workflows are combinations of rules that serve a specific purpose: 
quality check of sequencing data, ChIP-seq peaks analysis...


Workflows (.wf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File extension: *.wf

*todo*

Rules (.rules)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


*todo*

Python scripts (.py)
----------------------------------------------------------------

*todo*

R scripts
----------------------------------------------------------------

*todo*


Configuration files (yaml)
----------------------------------------------------------------

*todo*


R markdown files (.Rmd)
----------------------------------------------------------------

*todo*

Tabulated files (.tab)
----------------------------------------------------------------

We use tabulated files in order to define and describe the samples 
to be processed in the workflows. 

Examples of these files are available in the *examples* folder of the 
library. 

Sample description files (samples.tab)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*todo*

Experimental design files (design.tab)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*todo*
