"""Snakemake workflow for the analysis of RNA-seq data.

This script runs the following steps:

- merging the reads from multiple lanes to obtain a single read file per sample (fastq)
- read quality
- read mapping onto a reference genome
- counting the reads per gene

The parameters should be specified in a configuration file (yaml or json).

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow.py \
        --configfile [path_to_your_config_file.yaml] \
        [targets]

Flowcharts:
    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/Athaliana.py \
        --configfile [path_to_your_config_file.yaml] \
        --force flowcharts


Authors: Justine Long, Jeanne Cheneby, Lucie Khamvongsa, Claire Rioualen & Jacques van Helden
Contact: Jacques.van-Helden@univ-amu.fr
"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time#rm?
import datetime
import pandas as pd


#================================================================#
#                        Configuration
#================================================================#

## Config file must be specified on the command line, with the option --configfile

workdir: config["dir"]["base"]

# Beware: verbosity messages are incompatible with the flowcharts
verbosity = int(config["verbosity"])

# #================================================================#
# # Define suffixes for each step of the workflow. Note: this could
# # alternatively be done in the config file but we would then not be
# # able to build suffixes from other config values due to JSON
# # limitations.
# # ================================================================#
# config["suffix"]["trimmed"] = "sickle_pe_q" + config["sickle"]["threshold"]
# config["suffix"]["mapped"] = config["suffix"]["trimmed"] + "_bowtie2_pe"
# config["suffix"]["featurecounts"] = config["suffix"]["mapped"] + "_featurecounts"
# config["suffix"]["sorted_pos"] = config["suffix"]["mapped"] + "_sorted_pos"
# config["suffix"]["sorted_name"] = config["suffix"]["mapped"] + "_sorted_name"
# config["suffix"]["htseq_counts"] = config["suffix"]["sorted_name"] + "_HTSeqcount"
# config["suffix"]["deg"] = "sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"]
# config["suffix"]["edgeR"] = config["suffix"]["deg"] + config["edgeR"]["suffix"]
# config["suffix"]["DESeq2"] = config["suffix"]["deg"] + config["DESeq2"]["suffix"]
#
# #================================================================#
# # Define global variables
# #================================================================#
# NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
#
#================================================================#
# Import snakemake rules and python utilities
#================================================================#
if not ("dir" in config.keys()) & ("fg_lib" in config["dir"].keys()) :
    sys.exit("The parameter config['dir']['fg_lib'] should be specified in the config file.")

FG_LIB = os.path.abspath(config["dir"]["fg_lib"])
RULES = os.path.join(FG_LIB, "scripts/snakefiles/rules")
PYTHON = os.path.join(FG_LIB, "scripts/snakefiles/python_lib")

include: os.path.join(PYTHON, "util.py")                        ## Python utilities for our snakemake workflows
# include: os.path.join(RULES, "util.rules")                    ## Snakemake utilities
include: os.path.join(RULES, "merge_lanes.rules")               ## Build genome index for bowtie2 (read mapping with gaps)
include: os.path.join(RULES, "fastqc.rules")                    ## Quality control with fastqc
include: os.path.join(RULES, "count_reads.rules")               ## Count reads in different file formats
include: os.path.join(RULES, "bowtie2_build.rules")             ## Build genome index for bowtie2 (read mapping with gaps)
include: os.path.join(RULES, "bowtie2_paired_ends.rules")     ## Paired-ends read mapping with bowtie version 2 (support gaps)
# include: os.path.join(RULES, "sickle_paired_ends.rules")        ## Trimming with sickle
# include: os.path.join(RULES, "flowcharts.rules")              ## Draw flowcharts (dag and rule graph)
# include: os.path.join(RULES, "genome_coverage.rules")         ## Compute density profiles in bedgraph format
# include: os.path.join(RULES, "htseq.rules")                   ## Count reads per gene with htseq-count
# include: os.path.join(RULES, "featurecounts.rules")           ## Count reads per gene with R subread::featurecounts
#
#

#================================================================#
#                      Data & wildcards                          #
#================================================================#

# # Raw data
# READS = config["dir"]["reads_source"]

#----------------------------------------------------------------#
# Read sample descriptions
#----------------------------------------------------------------#

# Read the sample description file
SAMPLE_DESCR = read_table(config["files"]["samples"], verbosity=verbosity)
SAMPLE_IDS = SAMPLE_DESCR.iloc[:,0] ## First column MUST contain the sample ID
SAMPLE_CONDITIONS = SAMPLE_DESCR['condition'] ## Second column MUST contain condition for each sample
SAMPLE_NAMES = SAMPLE_DESCR['title'] ## Sample-wise label
SAMPLE_DIRS = SAMPLE_DESCR['folder']
FASTQ_R1 = SAMPLE_DESCR['fastq_R1']
FASTQ_R2 = SAMPLE_DESCR['fastq_R2']
FASTQ = list(FASTQ_R1) + list(FASTQ_R2)

# Verbosity
if (verbosity >= 1):
    print("Sample descriptions:\t" + config["files"]["samples"])
    if (verbosity >= 3):
        print("\tSample IDs:\t" + ";".join(SAMPLE_IDS))
        print("\tConditions:\t" + ";".join(SAMPLE_CONDITIONS))
        print("\tSample names:\t" + ";".join(SAMPLE_NAMES))
        print("\tSample folders:\t" + ";".join(SAMPLE_DIRS))


# #----------------------------------------------------------------#
# # Merge lanes per sample
# #----------------------------------------------------------------#
# rule merge_lanes:
#     """
#     Merge lanes (fastq) of the same sample and end in a single fastq file.
#
#     ince the file naming conventions are highly dependent on the sequencing
#     platform, the file grouping is read from a user-provided text file with
#     tab-separated values (extension .tsv). This file must have been specified
#     in the config file, as config["files"]["lane_merging"].
#
#     This file must contain at least two columns with this precise header:
#         source_file
#         merged_file
#
#     There should be a N to 1 correspondence from source file to merge file
#     (each source file should in principle be assigned to a single merged file).
#
#     Source files are supposed to be compressed fastq sequence files (.fastq.gz).
#
#     The output file is an uncompressed fastq file, because bowtie version 1
#     does not support gzipped files as input.
#
#     """
#     input: config["files"]["lane_merging"]
#     # output: config["dir"]["results"] + "_lane_merging_benchmark.json"
#     log: config["dir"]["results"] + "_lane_merging_log.txt"
#     benchmark: config["dir"]["results"] + "_lane_merging_benchmark.json"
#     run:
#         if (verbosity >= 1):
#             print("Lane merging table:\t" + config["files"]["lane_merging"])
#
#         # Read the lane merging table
#         lane_merging_table = read_table(config["files"]["lane_merging"], verbosity=verbosity)
#         source_file = lane_merging_table['source_file']
#         merged_file = lane_merging_table['merged_file']
#
#         # Build a dictionary indexed by merged file, where values are lists of files to be merged
#         merging_dict = {}
#         for s,m in zip(source_file, merged_file):
#             # print("\t".join([s,m]))
#             if (m in merging_dict):
#                 merging_dict[m].append(s)
#             else:
#                 merging_dict[m] = [s]
#
#         # Verbosity
#         if (verbosity >= 5):
#             print("\tsource_file:\t" + ";".join(source_file))
#             print("\tmerged_file:\t" + ";".join(merged_file))
#             print("\tmerging_dict:\t" + str(merging_dict))
#
#         # Merge the files
#         for m in merging_dict.keys():
#             # Check the output directory
#             m_dir = os.path.dirname(m)
#             if not os.path.exists(m_dir):
#                 os.makedirs(m_dir)
#
#             # Merge the source files
#             to_merge = merging_dict[m]
#             cmd = "gunzip -c " + " ".join(to_merge) + "> " + m
#             now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
#             if (verbosity >= 1):
#                 print(now + "\tMerging " + str(len(to_merge)) + " files into " + m)
#                 print("\t" + cmd)
#             os.system(cmd)


## Design
DESIGN = read_table(config["files"]["design"], verbosity=verbosity)
TREATMENT = DESIGN.iloc[:,0]
CONTROL = DESIGN.iloc[:,1]
if (verbosity >= 1):
    print("Design file:\t" + config["files"]["design"])
    if (verbosity >= 3):
        print("\tTREATMENT:\t" + ";".join(TREATMENT))
        print("\tCONTROL:\t" + ";".join(CONTROL))

# ## Ref genome
# GENOME = config["genome"]["version"]
#
# ## Results dir
# RESULTS_DIR = config["dir"]["results"]
# if not os.path.exists(RESULTS_DIR):
#     os.makedirs(RESULTS_DIR)

#================================================================#
# Define target file names
#================================================================#


#----------------------------------------------------------------#
# Genome index
#----------------------------------------------------------------#

GENOME_INDEX=config["bowtie2"]["index"] + "_benchmark.json"

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = [filename.replace('.fastq','_fastqc/') for filename in FASTQ]
RAW_READNB = [filename.replace('.fastq','_fastq_readnb.txt') for filename in FASTQ]

#----------------------------------------------------------------#
# Read mapping (alignment against reference genome)
#----------------------------------------------------------------#

## to avoid duplicates, fasta sequence should be moved to {genome} directly...
# BWA_INDEX = expand(config["dir"]["genomes"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
# BOWTIE2_INDEX = expand(config["dir"]["genomes"] + "{genome}/Bowtie2Index/{genome}.fa.1.bt2", genome=GENOME)
#
# MAPPING = expand(RESULTS_DIR + "{alignment}.sam", alignment=ALIGNMENT)
MAPPED_PE_SAM=expand(config["dir"]["sam"] + "/{sample_dir}/{sample_id}_bowtie2_pe.sam", zip, sample_dir=SAMPLE_IDS, sample_id=SAMPLE_IDS)
# MAPPED_PE_BAM=expand(config["dir"]["bam"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["mapped"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
if (verbosity >= 0):
    print("\tMAPPED_PE_SAM:\t" + ";".join(MAPPED_PE_SAM))


# #----------------------------------------------------------------#
# # Raw reads
# #----------------------------------------------------------------#
#
# # THIS SHOULD BE IMPROVED: for the time being I collect the file
# # names by listing the files corresponding to a given pattern. THE
# # INPUT FILE NAMES SHOULD BE PROVIDED IN THE SAMPLE DESCRIPTION FILE
# # !!!
# RAWR_L1R1, RAWR_L1R1_DIRS, RAWR_L1R1_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads_source"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
# RAWR_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
# RAWR_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + ".fastq", zip, sample_dir=RAWR_L1R1_DIRS, sample_basename=RAWR_L1R1_BASENAMES)
# RAWR_MERGED=RAWR_MERGED_FWD + RAWR_MERGED_REV
#
# # List separately the FORWARD and REVERSE raw read files
# RAWR_FILES_FWD, RAWR_DIRS_FWD, RAWR_BASENAMES_FWD=glob_multi_dir(SAMPLE_DIRS, "*" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_fwd"] + ".fastq.gz")
# RAWR_FILES_REV, RAWR_DIRS_REV, RAWR_BASENAMES_REV=glob_multi_dir(SAMPLE_DIRS, "*" + config["suffix"]["reads_rev"] + ".fastq.gz", config["dir"]["reads"], config["suffix"]["reads_rev"] + ".fastq.gz")
#
# # List all the raw read files
# RAWR_FILES, RAWR_DIRS, RAWR_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_R*_001.fastq.gz", config["dir"]["reads"], ".fastq.gz")
#
# #----------------------------------------------------------------#
# # Trimmed reads
# #----------------------------------------------------------------#
#
# # Merge trimmed reads. Note: I use a trick to obtain one directory
# # name per group of lanes: I only glob the first lane, and I use the
# # list of directories and basenames.
# SAMPLE_L1R1, PAIRED_DIRS, PAIRED_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["reads"], "_L001" + config["suffix"]["reads_fwd"] + ".fastq.gz")
# TRIMMED_MERGED_FWD=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_fwd"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# TRIMMED_MERGED_REV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged" + config["suffix"]["reads_rev"] + "_sickle_pe_q" + config["sickle"]["threshold"] + ".fastq", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# TRIMMED_MERGED=TRIMMED_MERGED_FWD + TRIMMED_MERGED_REV
# if verbosity >= 3:
#     print ("PAIRED_DIRS:\n\t" + "\n\t".join(PAIRED_DIRS))
#     print ("PAIRED_BASENAMES:\n\t" + "\n\t".join(PAIRED_BASENAMES))
#
# # Trimmed reads
# #TRIMMED_SUMMARIES = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_summary.txt", zip, reads=RAWR_BASENAMES_FWD, sample_dir=RAWR_DIRS_FWD)
# TRIMMED_FILES, TRIMMED_DIRS, TRIMMED_BASENAMES=glob_multi_dir(SAMPLE_DIRS, "*_R*_001_trimmed_thr" + config["sickle"]["threshold"] + ".fastq.gz", config["dir"]["reads"], ".fastq.gz")
#
# #----------------------------------------------------------------#
# # Quality control
# #----------------------------------------------------------------#
#
# #RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)
#
# MERGED_RAWR_QC_FWD = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
# MERGED_RAWR_QC_REV = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
# MERGED_RAWR_QC = MERGED_RAWR_QC_FWD + MERGED_RAWR_QC_REV
# #MERGED_RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
# #MERGED_RAWR_PREFIXES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
# #MERGED_RAWR_QC = expand("{fastq_prefix}_fastqc", fastq_prefix=MERGED_PREFIXES)
#
# TRIMMED_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=TRIMMED_BASENAMES, sample_dir=TRIMMED_DIRS)
# if verbosity >= 3:
#     print("TRIMMED_QC\t" + ";".join(TRIMMED_QC))
#
# #----------------------------------------------------------------#
# # Mapped reads
#
# #MAPPED_FILES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_bowtie2.sam", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)
#
# # Note: I use bowtie2 (paired-end with gaps) although we don't need
# # gaps in this project, because the sequence counts with htseq-count
# # was causing problems with bowtie1. I should revise this at some
# # point.
#
# MAPPED_PE_SAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["mapped"] + ".sam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# MAPPED_PE_BAM=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["mapped"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# MAPPED_PE_SORTED=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_pos"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# MAPPED_PE_SORTED_BY_NAME=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_name"] + ".bam", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# if (verbosity >= 3):
#     print ("MAPPED_PE_SORTED_BY_NAME:\n\t" + "\n\t".join(MAPPED_PE_SORTED_BY_NAME))
# GENOMECOV=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_pos"] + "_genomecov.tdf", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# GENOMECOV_PLUS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_pos"] + "_genomecov_strand+.tdf", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# GENOMECOV_MINUS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["sorted_pos"] + "_genomecov_strand-.tdf", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
#
# #----------------------------------------------------------------#
# # Read counts per gene (done with htseq-count)
# #----------------------------------------------------------------#
#
# # THERE SEEMS TO BE A PROBLEM WITH HTSEQ AND POSITION-SORTED BAM
# # FILES. This has been discussed on SesqAnswers:
# # http://seqanswers.com/forums/showthread.php?t=41531
#
# HTSEQ_COUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["htseq_counts"] + ".tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# if (verbosity >= 3):
#     print ("HTSEQ_COUNTS:\n\t" + "\n\t".join(HTSEQ_COUNTS))
#
# # Since the program Subreads featureCounts is MUCH faster (30 times)
# # than htseq-count, and does not required bam sorting, I switch to
# # featureCounts.
#
# FEATURECOUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["featurecounts"] + ".tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
# if (verbosity >= 3):
#     print ("FEATURECOUNTS:\n\t" + "\n\t".join(FEATURECOUNTS))
#
# #COUNT_FILES=HTSEQ_COUNTS
# COUNT_FILES=FEATURECOUNTS
#
# #----------------------------------------------------------------#
# # Differential expression analysis
# #----------------------------------------------------------------#
#
# PARAMS_R = config["dir"]["results"] + "/DEG/" + config["suffix"]["deg"] + "_params.R"
# ALL_COUNTS = config["dir"]["results"] + "/DEG/" + config["suffix"]["deg"] + "_all_counts.tab"
# if (verbosity >= 2):
#     print ("PARAMS_R:\t" + PARAMS_R)
#     print ("ALL_COUNTS:\t" + ALL_COUNTS)
#     print("COUNT_FILES\t" + ";".join(COUNT_FILES))
#
# #================================================================#
# # Rule definitions
# #================================================================#
#
# # Note: these rules must be loaded after having defined some global
# # variables COUNT_FILES, PARAMS_R, ALL_COUNTS.
# include: os.path.join(RULES, "allcount_params.rules")   ## Produce the count table from sample-based count files + the parameters for differential analysis
#
#
# # Read the analysis design file
# DESIGN = read_table(config["files"]["analyses"], verbosity=verbosity)
# config["Diff_Exp"]["cond1"] = COND_1 = DESIGN['cond1']
# config["Diff_Exp"]["cond2"] = COND_2 = DESIGN['cond2']
# if (verbosity >= 1):
#     print("Analysis design:\t" + config["files"]["analyses"])
#     if (verbosity >= 2):
#         print("\tCondition 1:\t" + ";".join(COND_1))
#         print("\tCondition 2:\t" + ";".join(COND_2))
#
# # Detect differentially expressed genes wit edgeR
# RESULTS_EDGER = expand(config["dir"]["results"] + "/DEG/{cond_1}_vs_{cond_2}/{cond_1}_vs_{cond_2}_" + config["suffix"]["edgeR"] +".tab", zip, cond_1=COND_1, cond_2=COND_2)
#
# include: os.path.join(RULES, "diff_expr.rules")                  ## Differential expression analysis with BioConductor edgeR and DESeq2 packates
# #include: os.path.join(RULES, "edgeR.rules")                  ## Differential expression analysis with BioConductor edgeR package
# #include: os.path.join(RULES, "DESeq2.rules")                 ## Differential expression analysis with BioConductor DESeq2 package
#
#
#================================================================#
#                        Rule all                                #
#================================================================#

rule all:
	"""
	Run all the required analyses.
	"""
	input: GENOME_INDEX, RAW_QC, RAW_READNB, MAPPED_PE_SAM
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

# #----------------------------------------------------------------#
# # Get all targets
# rule all:
#     """
#     Run all the required analyses
#     """
# #    input: TRIMMED_SUMMARIES ## Still working ?
# #    input: MERGED_RAWR_QC, RAWR_MERGED, TRIMMED_MERGED, TRIMMED_QC, MAPPED_PE_SAM, MAPPED_PE_BAM,
#     input: GENOMECOV, GENOMECOV_PLUS, GENOMECOV_MINUS
# #        MAPPED_PE_SORTED,  \
# #        GENOMECOV, \
# #        HTSEQ_COUNTS, \
# #        FEATURECOUNTS, \
# #        HTSEQ_COUNTS, \
# #        ALL_COUNTS, \
# #        RESULTS_EDGER
#     params: qsub=config["qsub"]
#     shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"
#
# ruleorder: sickle_paired_ends > merge_lanes
# ruleorder: bowtie2_paired_end > merge_lanes
#
# # ruleorder: bowtie_paired_end > merge_lanes

#
#
# #----------------------------------------------------------------#
# # Merge lanes per sample
# #----------------------------------------------------------------#
# rule merge_lanes:
#     """
#     Merge lanes of the same sample and end in a single file.  The input
#     files are compressed (.fastq.gz) but the output file is in
#     uncompressed fastq format, because bowtie version 1 does not
#     support gzipped files as input.
#
#     """
#     input: L1 = "{reads_prefix}_L001_{reads_suffix}.fastq.gz", \
#         L2 = "{reads_prefix}_L002_{reads_suffix}.fastq.gz", \
#         L3 = "{reads_prefix}_L003_{reads_suffix}.fastq.gz", \
#         L4 = "{reads_prefix}_L004_{reads_suffix}.fastq.gz"
#     output: "{reads_prefix}_merged_{reads_suffix}.fastq"
#     log: "{reads_prefix}_merged_{reads_suffix}.log"
#     benchmark: "{reads_prefix}_merged_{reads_suffix}_benchmark.json"
#     params: qsub = config["qsub"] + " -e {reads_prefix}__merged_{reads_suffix}_qsub.err  -o {reads_prefix}__merged_{reads_suffix}_qsub.out"
#     shell: "gunzip -c {input.L1} {input.L2} {input.L3} {input.L4} | gzip > {output}"

# #----------------------------------------------------------------#
# # Build the report (including DAG and rulegraph flowcharts).
# from snakemake.utils import report
#
# # Bulleted list of samples for the report
# SAMPLE_DIRS_OL=report_numbered_list(SAMPLE_DIRS)
# RAWR_MERGED_OL=report_numbered_list(RAWR_MERGED)
# TRIMMED_MERGED_OL=report_numbered_list(TRIMMED_MERGED)
# MAPPED_PE_SAM_OL=report_numbered_list(MAPPED_PE_SAM)
# MAPPED_PE_BAM_OL=report_numbered_list(MAPPED_PE_BAM)
# MAPPED_PE_SORTED_OL = report_numbered_list(MAPPED_PE_SORTED)
# MAPPED_PE_SORTED_BY_NAME_OL = report_numbered_list(MAPPED_PE_SORTED_BY_NAME)
# HTSEQ_COUNTS_OL = report_numbered_list(HTSEQ_COUNTS)
# FEATURECOUNTS_OL = report_numbered_list(FEATURECOUNTS)
# COUNT_FILES_OL = report_numbered_list(COUNT_FILES)
#
# rule report:
#     """
#     Generate a report with the list of datasets + summary of the results.
#     """
#     input:  dag=config["dir"]["reports"] + "/" + "dag.pdf", \
#             dag_png=config["dir"]["reports"] + "/" + "dag.png", \
#             rulegraph=config["dir"]["reports"] + "/" + "rulegraph.pdf", \
#             rulegraph_png=config["dir"]["reports"] + "/" + "rulegraph.png", \
#             all_counts = ALL_COUNTS, \
#             params_r = PARAMS_R
#     output: html=config["dir"]["reports"] + "/report.html"
#     run:
#         report("""
#         =================================
#         RNA-seq analysis - Béatrice Roche
#         =================================
#
#         :Date:                 {NOW}
#         :Project:              Béatrice Roche
#         :Analysis workflow:    Jacques van Helden
#
#         Contents
#         ========
#
#         - `Flowcharts`_
#         - `Datasets`_
#              - `Sample directories`_
#              - `Raw reads`_
#              - `Trimmed`_
#              - `Mapped`_
#              - `Count files`_
#
#         -----------------------------------------------------
#
#         Flowcharts
#         ==========
#
#         - Sample treatment: dag_
#         - Workflow: rulegraph_
#
#         .. image:: rulegraph.png
#
#         -----------------------------------------------------
#
#         Datasets
#         ========
#
#         Sample directories
#         ------------------
#
#         {SAMPLE_DIRS_OL}
#
#         Raw reads
#         ---------
#
#         (merged lanes per sample)
#
#         {RAWR_MERGED_OL}
#
#         Trimmed
#         -------
#
#         {TRIMMED_MERGED_OL}
#
#         Mapped
#         ------
#
#         Sam format (uncompressed)
#
#         {MAPPED_PE_SAM_OL}
#
#         Bam format (compressed)
#
#         {MAPPED_PE_BAM_OL}
#
#         Bam format (sorted by positions)
#
#         {MAPPED_PE_SORTED_OL}
#
#         Bam format (sorted by names)
#
#         {MAPPED_PE_SORTED_BY_NAME_OL}
#
#         Count files
#         -----------
#
#         htseq-count results (paired-ends, no multi overlap)
#
#         {HTSEQ_COUNTS_OL}
#
#         Subread featureCounts results (multi-overlaps)
#
#         ! temporarily: paired-ends option *inactive* due to problem
#
#         {FEATURECOUNTS_OL}
#
#         Count files for differential expression analysis
#
#         {COUNT_FILES_OL}
#
#         Count table
#         -----------
#
#         - Count table (one row per gene, one column per sample): all_counts_
#
#         R parameters
#         ------------
#
#         Parameters passed to R for differential expression anlysis
#
#         params_r_
#
#         -----------------------------------------------------
#
#         """, output.html, metadata="Jacques van Helden (Jacques.van-Helden@univ-amu.fr)", **input)
#
#
# # TO CHECK
# #   https://github.com/leipzig/snakemake-example/blob/master/Snakefile
# #   Report generated with R Sweave
