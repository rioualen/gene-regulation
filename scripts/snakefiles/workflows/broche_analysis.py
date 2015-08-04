"""Workflow for the analysis of transcriptomic RNA-seq data for
Beatrice Roche's project, treated by the TGML platform on June 2015.

Reference genome: Escherichia_coli_str_k_12_substr_mg1655

Usage: snakemake -p -c "qsub {params.qsub}" -j 30

"""

#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
import datetime
from snakemake.utils import R
configfile: "/home/jvanheld/fg-chip-seq/scripts/snakefiles/workflows/broche_analysis_config.json"
workdir: config["dir"]["base"]
#workdir: os.getcwd() # Local Root directoray for the project. Should be adapted for porting.

################################################################
## Define global variables
################################################################
NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

################################################################
## Import snakemake rules and python library
include: config["dir"]["fg-rules"] + "/util.py"                   ## Python utilities
include: config["dir"]["fg-rules"] + "/util.rules"                ## Snakemake utilities
include: config["dir"]["fg-rules"] + "/count_reads.rules"         ## Count reads in different file formats
include: config["dir"]["fg-rules"] + "/fastqc.rules"              ## Quality control with fastqc
include: config["dir"]["fg-rules"] + "/flowcharts.rules"          ## Draw flowcharts (dag and rule graph)
include: config["dir"]["fg-rules"] + "/sickle_paired_ends.rules"  ## Trimming with sickle
#include: config["dir"]["fg-rules"] + "/bowtie_build.rules"        ## Read mapping with bowtie version 1 (no gap)
#include: config["dir"]["fg-rules"] + "/bowtie_paired_ends.rules"  ## Paired-ends read mapping with bowtie version 1 (no gap)
include: config["dir"]["fg-rules"] + "/bowtie2_build.rules"        ## Read mapping with bowtie version 2 (suports gaps)
include: config["dir"]["fg-rules"] + "/bowtie2_paired_ends.rules"  ## Paired-ends read mapping with bowtie version 2 (support gaps)
include: config["dir"]["fg-rules"] + "/htseq.rules"        ## Count reads per gene with htseq-count
include: config["dir"]["fg-rules"] + "/featurecounts.rules"        ## Count reads per gene with R subread::featurecounts

################################################################
## Define the lists of requested files
################################################################

## Read the list of sample IDs from the sample description file
SAMPLE_IDS=read_sample_ids(config["files"]["samples"])

## List the merged raw reads
RAWR_L1R1, RAWR_L1R1_DIRS, RAWR_L1R1_BASENAMES=glob_multi_dir(SAMPLE_IDS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads_source"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
RAWR_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
RAWR_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
RAWR_MERGED=RAWR_MERGED_FWD + RAWR_MERGED_REV

## List separately the FORWARD and REVERSE raw read files, which will be submitted to quality control
RAWR_FILES_FWD, RAWR_DIRS_FWD, RAWR_BASENAMES_FWD=glob_multi_dir(SAMPLE_IDS, "*" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_fwd"] + ".fastq.gz")
RAWR_FILES_REV, RAWR_DIRS_REV, RAWR_BASENAMES_REV=glob_multi_dir(SAMPLE_IDS, "*" + config["suffix"]["reads_rev"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_rev"] + ".fastq.gz")


## Merge trimmed reads I use a trick to obtain one directory name per
## group of lanes: I only glob the first lane, and I use the list of
## directories and basenames
SAMPLE_L1R1, PAIRED_DIRS, PAIRED_BASENAMES=glob_multi_dir(SAMPLE_IDS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
TRIMMED_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
TRIMMED_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
TRIMMED_MERGED=TRIMMED_MERGED_FWD + TRIMMED_MERGED_REV

## Bowtie version 1: paired-end read mapping without gap
## Bowtie version 2: paired-end read mapping with gap
MAPPED_PE_SAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe.sam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_BAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe.bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_SORTED=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_pos.bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)

MAPPED_PE_SORTED_BY_NAME=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_name.bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)

## Problem with HTSeq and position-sorted bam files ? To be chacked later
## HTSEQ_COUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_pos_HTSeqcount.tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
HTSEQ_COUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_name_HTSeqcount.tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)

## List all the raw read files, which will be submitted to quality control
RAWR_FILES, RAWR_DIRS, RAWR_BASENAMES=glob_multi_dir(SAMPLE_IDS, "*_R*_001.fastq.gz", config["dir"]["reads"], ".fastq.gz")

## Quality control
#RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)
MERGED_RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
#MERGED_RAWR_PREFIXES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
#MERGED_RAWR_QC = expand("{fastq_prefix}_fastqc", fastq_prefix=MERGED_PREFIXES)

## Trimmed reads
# TRIMMED_SUMMARIES = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_summary.txt", zip, reads=RAWR_BASENAMES_FWD, sample_dir=RAWR_DIRS_FWD)
# TRIMMED_FILES, TRIMMED_DIRS, TRIMMED_BASENAMES=glob_multi_dir(SAMPLE_IDS, "*_R*_001_trimmed_thr" + config["sickle"]["threshold"] + ".fastq.gz", config["dir"]["reads"], ".fastq.gz")
# TRIMMED_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=TRIMMED_BASENAMES, sample_dir=TRIMMED_DIRS)

## Mapped reads
#MAPPED_FILES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_bowtie2.sam", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)

rule all: 
    """
    Run all the required analyses
    """
#    input: TRIMMED_SUMMARIES, TRIMMED_QC
#    input: MERGED_RAWR_QC, RAWR_MERGED, TRIMMED_MERGED, MAPPED_PE_SAM, MAPPED_PE_BAM, MAPPED_PE_SORTED
    input: MAPPED_PE_SAM, MAPPED_PE_BAM, MAPPED_PE_SORTED
#    input: MAPPED_PE_SORTED_BY_NAME, HTSEQ_COUNTS
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
SAMPLE_IDS_OL=report_numbered_list(SAMPLE_IDS)
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

        {SAMPLE_IDS_OL} 

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
