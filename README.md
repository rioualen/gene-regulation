# Shared code for France Genomique WP2.6 - Regulation

## Purpose of this repository

This git repository holds shared code for the analysis of Next
Generation Sequencing data related to gene regulation: ChIP-seq,
RNA-seq, and related technologies.

The goals are multiple.

1. Avoid duplication of efforts by sharing the developments of
different bioinformaticians involved in ChIP-seq and RNA-seq projects.

2. Ensure portability and re-usability of the code.

3. Enable validation of the code by independent users.


**Beware!** This folder should only contain re-usable code: 

- R libraries,
- snakemake rules, 
- generic snakemake workflows than can be included as sub-workflows in
  project-specific snakefiles.

Project-specific code should be handled in user-specific directories,
and should be no means be submitted to the distributed git repository.

## Contact

- Claire Rioualen <claire.rioualen@inserm.fr>
- Jacques van Helden <Jacques.van-helden@univ-amu.fr>

