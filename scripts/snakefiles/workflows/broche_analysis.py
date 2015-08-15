"""Workflow for the analysis of transcriptomic RNA-seq data for
Beatrice Roche's project, treated by the TGML platform on June 2015.

Reference genome: Escherichia_coli_str_k_12_substr_mg1655

Usage: snakemake -p -c "qsub {params.qsub}" -j 30

Analyse différentielle:

Il s'agit de deux projets RNAseq, A et B.  

- Pour le projet A, il s'agit de comparer chaque mutant (iscU, im et
  cyay) à la souche de référence (wt). Chaque souche a été faite en
  triplicat.

- Pour le projet B, il s'agit de comparer la souche 2 ( mutant IKE en
  triplicat sp5, 6, 7) à la souche 1 (mutant IK en triplicat sp1,2,3);
  puis la souche 3 (mutant stpA en triplicat sp8,9 10) à la souche 1;
  puis la souche 2 à la souche 3. Dans ce projet B il n'y a pas de
  souche wt car on veut comparer uniquement des mutants entre eux.

"""




#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
import datetime
from snakemake.utils import R
configfile: "scripts/snakefiles/workflows/broche_analysis_config.json"
workdir: config["dir"]["base"]
#workdir: os.getcwd() # Local Root directoray for the project. Should be adapted for porting.

## Beware: verbosity messages are incompatible with the flowcharts
verbosity = int(config["verbosity"])

################################################################
## Define global variables
################################################################
NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

################################################################
## Import snakemake rules and python library
include: config["dir"]["python_lib"] + "/util.py"                   ## Python utilities
include: config["dir"]["rules"] + "/util.rules"                     ## Snakemake utilities
include: config["dir"]["rules"] + "/count_reads.rules"              ## Count reads in different file formats
include: config["dir"]["rules"] + "/fastqc.rules"                   ## Quality control with fastqc
include: config["dir"]["rules"] + "/flowcharts.rules"               ## Draw flowcharts (dag and rule graph)
include: config["dir"]["rules"] + "/sickle_paired_ends.rules"       ## Trimming with sickle
#include: config["dir"]["rules"] + "/bowtie_build.rules"            ## Read mapping with bowtie version 1 (no gap)
#include: config["dir"]["rules"] + "/bowtie_paired_ends.rules"      ## Paired-ends read mapping with bowtie version 1 (no gap)
include: config["dir"]["rules"] + "/bowtie2_build.rules"            ## Read mapping with bowtie version 2 (suports gaps)
include: config["dir"]["rules"] + "/bowtie2_paired_ends.rules"      ## Paired-ends read mapping with bowtie version 2 (support gaps)
include: config["dir"]["rules"] + "/htseq.rules"                    ## Count reads per gene with htseq-count
include: config["dir"]["rules"] + "/featurecounts.rules"            ## Count reads per gene with R subread::featurecounts

################################################################
## Define the lists of requested files
################################################################

## Read the list of sample IDs from the sample description file
SAMPLE_IDS = read_sample_ids(config["files"]["samples"], column=1, verbosity=verbosity)
SAMPLE_TYPES = read_sample_ids(config["files"]["samples"], column=2, verbosity=verbosity)
SAMPLE_DIRS = read_sample_ids(config["files"]["samples"], column=3, verbosity=verbosity)

## Verbosity
if (verbosity >= 2):
    print ("Sample description file:\t" + config["files"]["samples"])
    print ("Sample IDs:\t" + ";".join(SAMPLE_IDS))
    print ("Sample directories:\t" + ";".join(SAMPLE_DIRS))

################################################################
## Raw read files.  

## THIS SHOULD BE IMPROVED: for the time being I collect the file
## names by listing the files corresponding to a given pattern. THE
## INPUT FILE NAMES SHOULD BE PROVIDED IN THE SAMPLE DESCRIPTION FILE
## !!!
RAWR_L1R1, RAWR_L1R1_DIRS, RAWR_L1R1_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads_source"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
RAWR_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
RAWR_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
RAWR_MERGED=RAWR_MERGED_FWD + RAWR_MERGED_REV

## List separately the FORWARD and REVERSE raw read files
RAWR_FILES_FWD, RAWR_DIRS_FWD, RAWR_BASENAMES_FWD=glob_multi_dir(SAMPLE_DIRS, "*" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_fwd"] + ".fastq.gz")
RAWR_FILES_REV, RAWR_DIRS_REV, RAWR_BASENAMES_REV=glob_multi_dir(SAMPLE_DIRS, "*" + config["suffix"]["reads_rev"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_rev"] + ".fastq.gz")

## List all the raw read files
RAWR_FILES, RAWR_DIRS, RAWR_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_R*_001.fastq.gz", config["dir"]["reads"], ".fastq.gz")



################################################################
## Trimmed reads

## Merge trimmed reads. Note: I use a trick to obtain one directory
## name per group of lanes: I only glob the first lane, and I use the
## list of directories and basenames.
SAMPLE_L1R1, PAIRED_DIRS, PAIRED_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
TRIMMED_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
TRIMMED_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
TRIMMED_MERGED=TRIMMED_MERGED_FWD + TRIMMED_MERGED_REV

## Trimmed reads
#TRIMMED_SUMMARIES = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_summary.txt", zip, reads=RAWR_BASENAMES_FWD, sample_dir=RAWR_DIRS_FWD)
TRIMMED_FILES, TRIMMED_DIRS, TRIMMED_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_R*_001_trimmed_thr" + config["sickle"]["threshold"] + ".fastq.gz", config["dir"]["reads"], ".fastq.gz")

################################################################
## Quality control

#RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)
MERGED_RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
#MERGED_RAWR_PREFIXES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
#MERGED_RAWR_QC = expand("{fastq_prefix}_fastqc", fastq_prefix=MERGED_PREFIXES)

TRIMMED_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=TRIMMED_BASENAMES, sample_dir=TRIMMED_DIRS)

################################################################
## Mapped reads

#MAPPED_FILES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_bowtie2.sam", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)

## Bowtie version 1: paired-end read mapping without gap.  
## Bowtie version 2: paired-end read mapping with gap.  
##
## Note: I finally opted for bowtie2, because the sequence counts wit
## htseq-count was causing problems with bowtie1. I should revise this
## at some point.

MAPPED_PE_SAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe.sam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_BAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe.bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_SORTED=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_pos.bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)

MAPPED_PE_SORTED_BY_NAME=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_name.bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
if (verbosity >= 3):
    print ("MAPPED_PE_SORTED_BY_NAME:\n\t" + "\n\t".join(MAPPED_PE_SORTED_BY_NAME))

################################################################
## Read counts per gene (done with htseq-count)

## THERE SEEMS TO BE A PROBLEM WITH HTSEQ AND POSITION-SORTED BAM
## FILES ?  TO BE CHACKED LATER

HTSEQ_COUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_name_HTSeqcount.tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
COUNT_FILES=HTSEQ_COUNTS
if (verbosity >= 3):
    print ("HTSEQ_COUNTS:\n\t" + "\n\t".join(HTSEQ_COUNTS))

################################################################
## Differential expression analysis
PARAMS_R = config["dir"]["results"] + "DEG/sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + "_params.R"
ALL_COUNTS = config["dir"]["results"] + "DEG/sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + "_allcounts.tab"
#OUT_HTSEQ = PARAMS_R, ALL_COUNTS
if (verbosity >= 0):
    print ("PARAMS_R:\t" + PARAMS_R)
    print ("ALL_COUNTS:\t" + ALL_COUNTS)
#    print ("OUT_HTSEQ:\t" + ",".join(OUT_HTSEQ))


################################################################
## Rules for differential analysis. 
##
## Note: these rules must be loaded after having defined some global
## variables COUNT_FILES, PARAMS_R, ALL_COUNTS.
include: config["dir"]["rules"] + "/HTseq_allcount_params.rules"    ## Produce the count table from sample-based count files + the parameters for differential analysis
#include: config["dir"]["rules"] + "/edgeR.rules"                   ## Differential expression analysis with BioConductor edgeR package
#include: config["dir"]["rules"] + "/DESeq2.rules"                  ## Differential expression analysis with BioConductor DESeq2 package


################################################################
## Targets
################################################################

rule all: 
    """
    Run all the required analyses
    """
#    input: TRIMMED_SUMMARIES, TRIMMED_QC
#    input: MAPPED_PE_SAM, MAPPED_PE_BAM, MAPPED_PE_SORTED
#    input: MAPPED_PE_SORTED_BY_NAME, HTSEQ_COUNTS
#    input: MERGED_RAWR_QC, RAWR_MERGED, TRIMMED_MERGED, TRIMMED_QC, MAPPED_PE_SAM, MAPPED_PE_BAM, MAPPED_PE_SORTED, MAPPED_PE_SORTED_BY_NAME, HTSEQ_COUNTS
    input: ALL_COUNTS
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

ruleorder: sickle_paired_ends > merge_lanes
ruleorder: bowtie2_paired_end > merge_lanes

# ruleorder: bowtie_paired_end > merge_lanes

rule merge_lanes:
    """
    Merge lanes of the same sample and end in a single file.  The input
    files are compressed (.fastq.gz) but the output file is in
    uncompressed fastq format, because bowtie version 1 does not
    support gzipped files as input.

    """
    input: L1 = "{reads_prefix}_L001_{reads_suffix}.fastq.gz", \
        L2 = "{reads_prefix}_L002_{reads_suffix}.fastq.gz", \
        L3 = "{reads_prefix}_L003_{reads_suffix}.fastq.gz", \
        L4 = "{reads_prefix}_L004_{reads_suffix}.fastq.gz"
    output: "{reads_prefix}_merged_{reads_suffix}.fastq" 
    log: "{reads_prefix}_merged_{reads_suffix}.log" 
    benchmark: "{reads_prefix}_merged_{reads_suffix}_benchmark.json" 
    params: qsub = config["qsub"] + " -e {reads_prefix}__merged_{reads_suffix}_qsub.err  -o {reads_prefix}__merged_{reads_suffix}_qsub.out"
    shell: "gunzip -c {input.L1} {input.L2} {input.L3} {input.L4} | gzip > {output}"

################################################################
## Build the report (including DAG and rulegraph flowcharts).
from snakemake.utils import report

## Bulleted list of samples for the report
SAMPLE_DIRS_OL=report_numbered_list(SAMPLE_DIRS)
RAWR_MERGED_OL=report_numbered_list(RAWR_MERGED)
TRIMMED_MERGED_OL=report_numbered_list(TRIMMED_MERGED)
MAPPED_PE_SAM_OL=report_numbered_list(MAPPED_PE_SAM)
MAPPED_PE_BAM_OL=report_numbered_list(MAPPED_PE_BAM)

rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  dag=config["dir"]["reports"] + "/" + "dag.pdf", \
            dag_png=config["dir"]["reports"] + "/" + "dag.png", \
            rulegraph=config["dir"]["reports"] + "/" + "rulegraph.pdf", \
            rulegraph_png=config["dir"]["reports"] + "/" + "rulegraph.png"
    output: html=config["dir"]["reports"] + "/report.html"
    run:
        report("""
        =================================
        RNA-seq analysis - Béatrice Roche
        =================================
        
        :Date:                 {NOW}
        :Project:              Béatrice Roche
        :Analysis workflow:    Jacques van Helden
        
        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets`_
             - `Sample directories`_
             - `Raw reads`_
             - `Trimmed`_
             - `Mapped`_

        -----------------------------------------------------

        Flowcharts
        ==========

        - Sample treatment: dag_
        - Workflow: rulegraph_

        .. image:: rulegraph.png

        -----------------------------------------------------

        Datasets
        ========
        
        Sample directories
        ------------------

        {SAMPLE_DIRS_OL} 

        Raw reads 
        ---------

        (merged lanes per sample)

        {RAWR_MERGED_OL}

        Trimmed
        -------

        {TRIMMED_MERGED_OL}

        Mapped
        ------

        Sam format (uncompressed)

        {MAPPED_PE_SAM_OL}

        Bam format (compressed)

        {MAPPED_PE_BAM_OL}

        -----------------------------------------------------

        """, output.html, metadata="Jacques van Helden (Jacques.van-Helden@univ-amu.fr)", **input)


## TO CHECK
##   https://github.com/leipzig/snakemake-example/blob/master/Snakefile
##   Report generated with R Sweave
