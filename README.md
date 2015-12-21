---
title: "France Genomique  Workpackage 2.6 - Gene Regulation
"
author: "Claire Rioualen"
date: "December 21, 2015"
output:
  html_document:
    fig_caption: yes
    toc: yes
    number_sections: true
  pdf_document:
    toc: yes
---

# Purpose of this repository

This git repository holds shared code for the analysis of Next
Generation Sequencing data related to gene regulation: ChIP-seq,
RNA-seq, and related technologies.

The goals are multiple.

1. Avoid duplication of efforts by sharing the developments of
different bioinformaticians involved in ChIP-seq and RNA-seq projects.

2. Ensure portability and re-usability of the code.

3. Enable validation of the code by independent users.

# Workflow examples

This repository contains basic workflow examples. 

Pre-requisites:

* R 3+
* Python 3.3+
* Snakemake 3.3+

See `doc/install_protocols/install_snakemake_workflows.Rmd` for details.

## Transcription factors

![alt text][factor]

### Execution 

The file `factor_workflow.py` is designed to perform ChIP-seq analysis on transcription factor studies. 

Two associated config files are designed to work with this workflow: `Athaliana-Myb.yml` and `Scerevisiae-GCN4.yml`.

In order to launch the workflows, you can enter the following commands from the root of the repository:

```
snakemake -s scripts/snakefiles/workflows/factor_workflow.py --configfile examples/Scerevisiae-GCN4/Scerevisiae-GCN4.yml -p -n
snakemake -s scripts/snakefiles/workflows/factor_workflow.py --configfile examples/Athaliana-Myb/Athaliana-Myb.yml -p -n
```
This ensures the workflows are executable. Remove the -n option to actually run them. 

### Dependencies

* [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
* [Bowtie 2](http://bowtie-bio.sourceforge.net/)
* [SAMtools](http://samtools.sourceforge.net/)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [bedtools](http://bedtools.readthedocs.org/)
* [HOMER](http://homer.salk.edu/homer/index.html)
* [MACS1.4](http://liulab.dfci.harvard.edu/MACS/index.html)
* [MACS2](https://github.com/taoliu/MACS/)
* [SWEMBL](http://www.ebi.ac.uk/~swilder/SWEMBL/)
* [RSAT](http://rsat.eu/)


## Histone marks

![alt text][histone]

This workflow is adapted for histone marks that show larger profiles than usual TF marks. 

### Execution 

The file `histone_workflow.py` is designed to perform ChIP-seq analysis on histone marks data. 

Two associated config files are designed to work with this workflow: `Celegans-H3K4me3.yml` and `Celegans-H3K27me3.yml`.

In order to launch the workflows, you can enter the following commands from the root of the repository:

```
snakemake -s scripts/snakefiles/workflows/histone_workflow.py --configfile examples/Celegans-H3K4me3/Celegans-H3K4me3.yml -p -n
snakemake -s scripts/snakefiles/workflows/histone_workflow.py --configfile examples/Celegans-H3K27me3/Celegans-H3K27me3.yml -p -n
```
This ensures the workflows are executable. Remove the -n option to actually run them. 

### Dependencies

* [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
* [BWA](http://bio-bwa.sourceforge.net/)
* [SAMtools](http://samtools.sourceforge.net/)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [bedtools](http://bedtools.readthedocs.org/)
* [HOMER](http://homer.salk.edu/homer/index.html)
* [MACS1.4](http://liulab.dfci.harvard.edu/MACS/index.html)
* [MACS2](https://github.com/taoliu/MACS/)
* [SWEMBL](http://www.ebi.ac.uk/~swilder/SWEMBL/)

## RNA-seq (TODO)

# Data organisation

From the example config files, it is assumed that, in order to run the workflows, you have the following file organisation:

```
/data/
    Athaliana-Myb/
        GSM1482283/
            SRR1554463.sra
        GSM1482284/
            SRR1554464.sra
        design.tab
        samples.tab
    Celegans-H3K27me3/
        GSM1217457/
            SRR1198569.sra
        GSM121745/
            SRR1198570.sra
        GSM121745/
            SRR1198571.sra
        GSM1217460/
            SRR1198572.sra
        design.tab
        samples.tab
    Celegans-H3K4me3/
        GSM120836/
            SRR952381.sra
        GSM1217457/
            SRR1198569.sra
        design.tab
        samples.tab
    Scerevisiae-GCN4/
        GSM147015/
            SRR1542426.sra
        GSM1470160/
            SRR1542427.sra
        GSM147016/
            SRR1542428.sra
        GSM147016/
            SRR1542429.sra
        GSM1470163/
            SRR1542430.sra
        GSM1470164/
            SRR1542431.sra
        design.tab
        samples.tab
    genomes/
        ce10/
            ce10.fa
        sacCer3/
            sacCer3.fa
        TAIR10/
            TAIR10.fa
```

All the data is available from the GEO platform. `samples.tab` and `design.tab` can be found in the repository examples section. 

# Documentation

More documentation can be found in the `doc` directory. (under construction)

It includes: 

* Instructions for the installation of dependencies (to be refhreshed)
* RSAT install guide (to be refreshed)
* Instructions for building a virtual machine on the IFB cloud or under VirtualBox
* Some wiki material



# Contact

- Claire Rioualen <claire.rioualen@inserm.fr>
- Jacques van Helden <Jacques.van-helden@univ-amu.fr>


[factor]: https://github.com/rioualen/gene-regulation/blob/master/examples/factor.png
[histone]: https://github.com/rioualen/gene-regulation/blob/master/examples/histone.png
