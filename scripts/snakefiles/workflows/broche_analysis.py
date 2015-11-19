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

DO TO
=====

- check why featurecounts takes ages with one particular sample when
  the option paired ends is active

- suppress bam sorted by positions

- compare paired end counts between htseq-counts and featurecounts,
  with the appropriate options (paired-ends, no multi overlap).



"""

#================================================================#
#                        Required Python libraries               #
#================================================================#

import os
import datetime
from snakemake.utils import R

#================================================================#
#                        Configuration
#================================================================#

configfile: "scripts/snakefiles/workflows/broche_analysis_config.json"
workdir: config["dir"]["base"]
#workdir: os.getcwd() # Local Root directoray for the project. Should be adapted for porting.

# Beware: verbosity messages are incompatible with the flowcharts
verbosity = int(config["verbosity"])

#================================================================#
# Define suffixes for each step of the workflow. Note: this could
# alternatively be done in the config file but we would then not be
# able to build suffixes from other config values due to JSON
# limitations.
# ================================================================#
config["suffix"]["trimmed"] = "sickle_pe_q" + config["sickle"]["threshold"]
config["suffix"]["mapped"] = config["suffix"]["trimmed"] + "_bowtie2_pe"
config["suffix"]["featurecounts"] = config["suffix"]["mapped"] + "_featurecounts"
config["suffix"]["sorted_pos"] = config["suffix"]["mapped"] + "_sorted_pos"
config["suffix"]["sorted_name"] = config["suffix"]["mapped"] + "_sorted_name"
config["suffix"]["htseq_counts"] = config["suffix"]["sorted_name"] + "_HTSeqcount"
config["suffix"]["deg"] = "sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"]
config["suffix"]["edgeR"] = config["suffix"]["deg"] + config["edgeR"]["suffix"]
config["suffix"]["DESeq2"] = config["suffix"]["deg"] + config["DESeq2"]["suffix"]

#================================================================#
# Define global variables
#================================================================#
NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

#================================================================#
# Import snakemake rules and python utilities
#================================================================#

include: config["dir"]["python_lib"] + "/util.py"                   ## Python utilities for our snakemake workflows
include: config["dir"]["rules"] + "/util.rules"                     ## Snakemake utilities
include: config["dir"]["rules"] + "/count_reads.rules"              ## Count reads in different file formats
include: config["dir"]["rules"] + "/fastqc.rules"                   ## Quality control with fastqc
include: config["dir"]["rules"] + "/flowcharts.rules"               ## Draw flowcharts (dag and rule graph)
include: config["dir"]["rules"] + "/sickle_paired_ends.rules"       ## Trimming with sickle
#include: config["dir"]["rules"] + "/bowtie_build.rules"            ## Read mapping with bowtie version 1 (no gap)
#include: config["dir"]["rules"] + "/bowtie_paired_ends.rules"      ## Paired-ends read mapping with bowtie version 1 (no gap)
include: config["dir"]["rules"] + "/bowtie2_build.rules"            ## Read mapping with bowtie version 2 (suports gaps)
include: config["dir"]["rules"] + "/bowtie2_paired_ends.rules"      ## Paired-ends read mapping with bowtie version 2 (support gaps)
include: config["dir"]["rules"] + "/genome_coverage.rules"          ## Compute density profiles in bedgraph format
include: config["dir"]["rules"] + "/htseq.rules"                    ## Count reads per gene with htseq-count
include: config["dir"]["rules"] + "/featurecounts.rules"            ## Count reads per gene with R subread::featurecounts


#================================================================#
# Read sample descriptions
#================================================================#

# Read the sample description file
SAMPLE_DESCR = read_table(config["files"]["samples"], verbosity=verbosity)
SAMPLE_DIRS = SAMPLE_DESCR['folder']
SAMPLE_IDS = SAMPLE_DESCR.iloc[:,0] ## First column MUST contain the sample ID
SAMPLE_CONDITIONS = SAMPLE_DESCR.iloc[:,1] ## Second column MUST contain condition for each sample

# Verbosity
if (verbosity >= 1):
    print("Sample descriptions:\t" + config["files"]["samples"])
    if (verbosity >= 2):
        print("\tSample IDs:\t" + ";".join(SAMPLE_IDS))
        print("\tSample folders:\t" + ";".join(SAMPLE_DIRS))
        print("\tConditions:\t" + ";".join(SAMPLE_CONDITIONS))

#================================================================#
# Define target file names
#================================================================#

#----------------------------------------------------------------#
# Raw reads
#----------------------------------------------------------------#

# THIS SHOULD BE IMPROVED: for the time being I collect the file
# names by listing the files corresponding to a given pattern. THE
# INPUT FILE NAMES SHOULD BE PROVIDED IN THE SAMPLE DESCRIPTION FILE
# !!!
RAWR_L1R1, RAWR_L1R1_DIRS, RAWR_L1R1_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads_source"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
RAWR_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
RAWR_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
RAWR_MERGED=RAWR_MERGED_FWD + RAWR_MERGED_REV

# List separately the FORWARD and REVERSE raw read files
RAWR_FILES_FWD, RAWR_DIRS_FWD, RAWR_BASENAMES_FWD=glob_multi_dir(SAMPLE_DIRS, "*" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_fwd"] + ".fastq.gz")
RAWR_FILES_REV, RAWR_DIRS_REV, RAWR_BASENAMES_REV=glob_multi_dir(SAMPLE_DIRS, "*" + config["suffix"]["reads_rev"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_rev"] + ".fastq.gz")

# List all the raw read files
RAWR_FILES, RAWR_DIRS, RAWR_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_R*_001.fastq.gz", config["dir"]["reads"], ".fastq.gz")

#----------------------------------------------------------------#
# Trimmed reads
#----------------------------------------------------------------#

# Merge trimmed reads. Note: I use a trick to obtain one directory
# name per group of lanes: I only glob the first lane, and I use the
# list of directories and basenames.
SAMPLE_L1R1, PAIRED_DIRS, PAIRED_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
TRIMMED_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
TRIMMED_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
TRIMMED_MERGED=TRIMMED_MERGED_FWD + TRIMMED_MERGED_REV
if verbosity >= 3:
    print ("PAIRED_DIRS:\n\t" + "\n\t".join(PAIRED_DIRS))
    print ("PAIRED_BASENAMES:\n\t" + "\n\t".join(PAIRED_BASENAMES))

# Trimmed reads
#TRIMMED_SUMMARIES = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_summary.txt", zip, reads=RAWR_BASENAMES_FWD, sample_dir=RAWR_DIRS_FWD)
TRIMMED_FILES, TRIMMED_DIRS, TRIMMED_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_R*_001_trimmed_thr" + config["sickle"]["threshold"] + ".fastq.gz", config["dir"]["reads"], ".fastq.gz")

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

#RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)

MERGED_RAWR_QC_FWD = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
MERGED_RAWR_QC_REV = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
MERGED_RAWR_QC = MERGED_RAWR_QC_FWD + MERGED_RAWR_QC_REV
#MERGED_RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
#MERGED_RAWR_PREFIXES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
#MERGED_RAWR_QC = expand("{fastq_prefix}_fastqc", fastq_prefix=MERGED_PREFIXES)

TRIMMED_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=TRIMMED_BASENAMES, sample_dir=TRIMMED_DIRS)
if verbosity >= 3:
    print("TRIMMED_QC\t" + ";".join(TRIMMED_QC))

#----------------------------------------------------------------#
# Mapped reads

#MAPPED_FILES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_bowtie2.sam", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)

# Note: I use bowtie2 (paired-end with gaps) although we don't need
# gaps in this project, because the sequence counts with htseq-count
# was causing problems with bowtie1. I should revise this at some
# point.

MAPPED_PE_SAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["mapped"] + ".sam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_BAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["mapped"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_SORTED=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_pos"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
MAPPED_PE_SORTED_BY_NAME=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_name"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
if (verbosity >= 3):
    print ("MAPPED_PE_SORTED_BY_NAME:\n\t" + "\n\t".join(MAPPED_PE_SORTED_BY_NAME))
GENOMECOV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_pos"] + "_genomecov.bedgraph", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)

#----------------------------------------------------------------#
# Read counts per gene (done with htseq-count)
#----------------------------------------------------------------#

# THERE SEEMS TO BE A PROBLEM WITH HTSEQ AND POSITION-SORTED BAM
# FILES. This has been discussed on SesqAnswers: 
# http://seqanswers.com/forums/showthread.php?t=41531

HTSEQ_COUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["htseq_counts"] + ".tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
if (verbosity >= 3):
    print ("HTSEQ_COUNTS:\n\t" + "\n\t".join(HTSEQ_COUNTS))

# Since the program Subreads featureCounts is MUCH faster (30 times)
# than htseq-count, and does not required bam sorting, I switch to
# featureCounts.

FEATURECOUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["featurecounts"] + ".tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
if (verbosity >= 3):
    print ("FEATURECOUNTS:\n\t" + "\n\t".join(FEATURECOUNTS))

#COUNT_FILES=HTSEQ_COUNTS
COUNT_FILES=FEATURECOUNTS

#----------------------------------------------------------------#
# Differential expression analysis
#----------------------------------------------------------------#

PARAMS_R = config["dir"]["results"] + "/DEG/" + config["suffix"]["deg"] + "_params.R"
ALL_COUNTS = config["dir"]["results"] + "/DEG/" + config["suffix"]["deg"] + "_all_counts.tab"
if (verbosity >= 2):
    print ("PARAMS_R:\t" + PARAMS_R)
    print ("ALL_COUNTS:\t" + ALL_COUNTS)
    print("COUNT_FILES\t" + ";".join(COUNT_FILES))

#================================================================#
# Rule definitions
#================================================================#

# Note: these rules must be loaded after having defined some global
# variables COUNT_FILES, PARAMS_R, ALL_COUNTS.
include: config["dir"]["rules"] + "/allcount_params.rules"    ## Produce the count table from sample-based count files + the parameters for differential analysis


# Read the analysis design file
DESIGN = read_table(config["files"]["analyses"], verbosity=verbosity)
config["Diff_Exp"]["cond1"] = COND_1 = DESIGN['cond1']
config["Diff_Exp"]["cond2"] = COND_2 = DESIGN['cond2']
if (verbosity >= 1):
    print("Analysis design:\t" + config["files"]["analyses"])
    if (verbosity >= 2):
        print("\tCondition 1:\t" + ";".join(COND_1))
        print("\tCondition 2:\t" + ";".join(COND_2))

# Detect differentially expressed genes wit edgeR
RESULTS_EDGER = expand(config["dir"]["results"] + "/DEG/{cond_1}_vs_{cond_2}/{cond_1}_vs_{cond_2}_" + config["suffix"]["edgeR"] +".tab", zip, cond_1=COND_1, cond_2=COND_2)

include: config["dir"]["rules"] + "/diff_expr.rules"                   ## Differential expression analysis with BioConductor edgeR and DESeq2 packates
#include: config["dir"]["rules"] + "/edgeR.rules"                   ## Differential expression analysis with BioConductor edgeR package
#include: config["dir"]["rules"] + "/DESeq2.rules"                  ## Differential expression analysis with BioConductor DESeq2 package


#----------------------------------------------------------------#
# Get all targets
rule all: 
    """
    Run all the required analyses
    """
#    input: TRIMMED_SUMMARIES ## Still working ?
#    input: MERGED_RAWR_QC, RAWR_MERGED, TRIMMED_MERGED, TRIMMED_QC, MAPPED_PE_SAM, MAPPED_PE_BAM, 
    input: MAPPED_PE_SORTED, GENOMECOV, HTSEQ_COUNTS, FEATURECOUNTS, ALL_COUNTS
    #input: HTSEQ_COUNTS, RESULTS_EDGER
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

ruleorder: sickle_paired_ends > merge_lanes
ruleorder: bowtie2_paired_end > merge_lanes

# ruleorder: bowtie_paired_end > merge_lanes


#----------------------------------------------------------------#
# Merge lanes per sample
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

#----------------------------------------------------------------#
# Build the report (including DAG and rulegraph flowcharts).
from snakemake.utils import report

# Bulleted list of samples for the report
SAMPLE_DIRS_OL=report_numbered_list(SAMPLE_DIRS)
RAWR_MERGED_OL=report_numbered_list(RAWR_MERGED)
TRIMMED_MERGED_OL=report_numbered_list(TRIMMED_MERGED)
MAPPED_PE_SAM_OL=report_numbered_list(MAPPED_PE_SAM)
MAPPED_PE_BAM_OL=report_numbered_list(MAPPED_PE_BAM)
MAPPED_PE_SORTED_OL = report_numbered_list(MAPPED_PE_SORTED)
MAPPED_PE_SORTED_BY_NAME_OL = report_numbered_list(MAPPED_PE_SORTED_BY_NAME)
HTSEQ_COUNTS_OL = report_numbered_list(HTSEQ_COUNTS)
FEATURECOUNTS_OL = report_numbered_list(FEATURECOUNTS)
COUNT_FILES_OL = report_numbered_list(COUNT_FILES)

rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  dag=config["dir"]["reports"] + "/" + "dag.pdf", \
            dag_png=config["dir"]["reports"] + "/" + "dag.png", \
            rulegraph=config["dir"]["reports"] + "/" + "rulegraph.pdf", \
            rulegraph_png=config["dir"]["reports"] + "/" + "rulegraph.png", \
            all_counts = ALL_COUNTS, \
            params_r = PARAMS_R
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
             - `Count files`_

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

        Bam format (sorted by positions)

        {MAPPED_PE_SORTED_OL}

        Bam format (sorted by names)

        {MAPPED_PE_SORTED_BY_NAME_OL}

        Count files
        -----------

        htseq-count results (paired-ends, no multi overlap)

        {HTSEQ_COUNTS_OL}

        Subread featureCounts results (multi-overlaps)

        ! temporarily: paired-ends option *inactive* due to problem

        {FEATURECOUNTS_OL}

        Count files for differential expression analysis

        {COUNT_FILES_OL}

        Count table
        -----------

        - Count table (one row per gene, one column per sample): all_counts_

        R parameters
        ------------

        Parameters passed to R for differential expression anlysis

        params_r_

        -----------------------------------------------------

        """, output.html, metadata="Jacques van Helden (Jacques.van-Helden@univ-amu.fr)", **input)


# TO CHECK
#   https://github.com/leipzig/snakemake-example/blob/master/Snakefile
#   Report generated with R Sweave
