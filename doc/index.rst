.. GeneRegulation documentation master file, created by
   sphinx-quickstart on Tue Mar 14 14:04:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
   
Gene Regulation
================================================================

This git repository holds shared code for the analysis of Next
Generation Sequencing data related to gene regulation: ChIP-seq,
RNA-seq, and related technologies.

It was developed as part of the France Genomique Workpackage 2.6: 
Gene expression regulation. 

One of the key goals is to ensure portability and re-usability of the code.

We have chosen to use the `Snakemake workflow management
system <https://bitbucket.org/snakemake/snakemake/wiki/Home>`__\ [1] in
order to build reproducible and flexible NGS analysis pipelines.


Contact
----------------------------------------------------------------

-  Claire Rioualen claire.rioualen@inserm.fr
-  Jacques van Helden Jacques.van-helden@univ-amu.fr

References
----------------------------------------------------------------

1. Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable
   bioinformatics workflow engine". Bioinformatics 2012.



User guide and reference
----------------------------------------------------------------


.. toctree::
    :numbered:
    :maxdepth: 3

    getting_started.rst
    gene-regulation_library.rst
    dependencies.rst
    tutorials.rst
    environments.rst
    snakemake.rst
    wiki.rst

..    other_tools.rst

Indices and tables
----------------------------------------------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


