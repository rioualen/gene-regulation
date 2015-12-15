# France Genomique  Workpackage 2.6 - Gene Regulation

Last update: 15/12/15

## Purpose of this repository

This git repository holds shared code for the analysis of Next
Generation Sequencing data related to gene regulation: ChIP-seq,
RNA-seq, and related technologies.

The goals are multiple.

1. Avoid duplication of efforts by sharing the developments of
different bioinformaticians involved in ChIP-seq and RNA-seq projects.

2. Ensure portability and re-usability of the code.

3. Enable validation of the code by independent users.

## Workflow examples

This repository contains basic workflow examples. 

### Transcription factors

The file `factor_workflow.py` is designed to perform ChIP-seq analysis on transcription factor studies. 

Two associated config files are designed to work with this workflow: `Athaliana-Myb.yml` and `Scerevisiae-GCN4.yml`.

In order to launch the workflows, you can enter the following commands from the root of the repository:

```
snakemake -s scripts/snakefiles/workflows/factor_workflow.py --configfile examples/Scerevisiae-GCN4/Scerevisiae-GCN4.yml -p -n
snakemake -s scripts/snakefiles/workflows/factor_workflow.py --configfile examples/Athaliana-Myb/Athaliana-Myb.yml -p -n
```
This ensures the workflows are executable. Remove the -n option to actually run them. 

### Histone marks

The file `histone_workflow.py` is designed to perform ChIP-seq analysis on histone marks data. 

Two associated config files are designed to work with this workflow: `Celegans-H3K4me3.yml` and `Celegans-H3K27me3.yml`.

In order to launch the workflows, you can enter the following commands from the root of the repository:

```
snakemake -s scripts/snakefiles/workflows/histone_workflow.py --configfile examples/Celegans-H3K4me3/Celegans-H3K4me3.yml -p -n
snakemake -s scripts/snakefiles/workflows/histone_workflow.py --configfile examples/Celegans-H3K27me3/Celegans-H3K27me3.yml -p -n
```
This ensures the workflows are executable. Remove the -n option to actually run them. 

### Data organisation

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

## Contact

- Claire Rioualen <claire.rioualen@inserm.fr>
- Jacques van Helden <Jacques.van-helden@univ-amu.fr>


