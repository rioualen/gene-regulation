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
include: os.path.join(RULES, "util.rules")                    ## Snakemake utilities
include: os.path.join(RULES, "flowcharts.rules")              ## Draw flowcharts (dag and rule graph)
include: os.path.join(RULES, "merge_lanes.rules")               ## Merge lanes by sample, based on a tab-delimited file indicating how to merge
include: os.path.join(RULES, "fastqc.rules")                    ## Quality control with fastqc
# include: os.path.join(RULES, "sickle_paired_ends.rules")        ## Trimming with sickle
include: os.path.join(RULES, "count_reads.rules")               ## Count reads in different file formats
#include: os.path.join(RULES, "bowtie2_build.rules")             ## Build genome index for bowtie2 (read mapping with gaps)
#include: os.path.join(RULES, "bowtie2_paired_ends.rules")       ## Paired-ends read mapping with bowtie version 2 (support gaps)
include: os.path.join(RULES, "subread_mapping_JvH.rules")       ## Read mapping with subreads
include: os.path.join(RULES, "genome_coverage.rules")         ## Compute density profiles in bedgraph format
# include: os.path.join(RULES, "htseq.rules")                   ## Count reads per gene with htseq-count
include: os.path.join(RULES, "featurecounts.rules")           ## Count reads per gene with R subread::featurecounts
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

# GENOME_INDEX=config["bowtie2"]["index"] + "_benchmark.json"
GENOME_INDEX=config["subread"]["index"] + ".files"

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = [filename.replace('.fastq','_fastqc/') for filename in FASTQ]
RAW_READNB = [filename.replace('.fastq','_fastq_readnb.txt') for filename in FASTQ]

#----------------------------------------------------------------#
# Read mapping (alignment against reference genome)
#----------------------------------------------------------------#

# Note: after having implemented rules for read-mapping with BWA,
# bowtie and bowtie2, we opted for subread-align, which is *much*
# faster than any of these.

## Aligned reads produced by subread-align (10 times faster than bowtie2).
SUBREADALIGN_PE_BAM=expand(config["dir"]["mapped_reads"] + "/{sample_dir}/{sample_id}_subread-align_pe.bam", zip, sample_dir=SAMPLE_IDS, sample_id=SAMPLE_IDS)
if (verbosity >= 3):
    print("\tSUBREADALIGN_PE_BAM:\t" + ";".join(SUBREADALIGN_PE_BAM))

## Genome coverage file (number of reads per genomic window), useful
## for visualisation. The bedgraph format is essentially used as
## intermediate to obtain TDF files, preferred by IGV.
SUBREADALIGN_PE_BG=expand(config["dir"]["mapped_reads"] + "/{sample_dir}/{sample_id}_subread-align_pe_sorted_pos.bg", zip, sample_dir=SAMPLE_IDS, sample_id=SAMPLE_IDS)
if (verbosity >= 3):
    print("\tSUBREADALIGN_PE_BG:\t" + ";".join(SUBREADALIGN_PE_BG))

## Genome coverage file (number of reads per genomic window), useful
## for visualisation. The TDF format is the standard for IGV.
SUBREADALIGN_PE_TDF=expand(config["dir"]["mapped_reads"] + "/{sample_dir}/{sample_id}_subread-align_pe_sorted_pos.tdf", zip, sample_dir=SAMPLE_IDS, sample_id=SAMPLE_IDS)
if (verbosity >= 3):
    print("\tSUBREADALIGN_PE_TDF:\t" + ";".join(SUBREADALIGN_PE_TDF))

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

#----------------------------------------------------------------#
# Read counts per gene (done with htseq-count)
#----------------------------------------------------------------#

# Since the program featureCounts (subread suite) is MUCH faster (30
# times) than htseq-count, and does not required bam sorting, I switch
# to featureCounts.


#FEATURECOUNTS=expand(config["dir"]["reads"] + "/{sample_dir}/{sample_basename}_merged_" + config["suffix"]["featurecounts"] + ".tab", zip, sample_dir=PAIRED_DIRS, sample_basename=PAIRED_BASENAMES)
COUNT_FILES=expand(config["dir"]["mapped_reads"] + "/{sample_dir}/{sample_id}_subread-align_pe_featurecounts.tab", zip, sample_dir=SAMPLE_IDS, sample_id=SAMPLE_IDS)
if (verbosity >= 3):
    print ("COUNT_FILES:\n\t" + "\n\t".join(COUNT_FILES))


## Print a tab-delimited file with the paths of count files per sample
if not ("files" in config.keys()) & ("count_file_paths" in config["files"].keys()) :
    sys.exit("The parameter config['files']['count_file_paths'] should be specified in the config file.")
## Check the directory for count file paths and create it if required
count_paths_dir = os.path.dirname(config['files']['count_file_paths'])
if not os.path.exists(count_paths_dir):
        os.makedirs(count_paths_dir)

COUNT_FILE_PATHS=pd.DataFrame({"SampleID": SAMPLE_IDS, "FilePath":COUNT_FILES})
COUNT_FILE_PATHS[["SampleID", "FilePath"]].to_csv(config["files"]["count_file_paths"], sep='\t', 
                        encoding='utf-8', header=True, index=False)
if (verbosity >= 1):
    print("Count file paths\t" + config["files"]["count_file_paths"])
    if (verbosity >= 3):
        print(COUNT_FILE_PATHS)

## Summary count table, with one row per gene and one column per sample
if not ("files" in config.keys()) & ("count_table" in config["files"].keys()) :
    sys.exit("The parameter config['files']['count_table'] should be specified in the config file.")
## Check the directory for count file paths and create it if required
count_table_dir = os.path.dirname(config['files']['count_table'])
if not os.path.exists(count_table_dir):
        os.makedirs(count_table_dir)
COUNT_TABLE=config['files']['count_table']


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
	input: GENOME_INDEX, RAW_QC, RAW_READNB, SUBREADALIGN_PE_BAM, SUBREADALIGN_PE_TDF, COUNT_FILES #, COUNT_TABLE
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

# #----------------------------------------------------------------#
# # Get all targets
# rule all:
#     """
#     Run all the required analyses
#     """
# #    input: TRIMMED_SUMMARIES ## Still working ?
# #    input: MERGED_RAWR_QC, RAWR_MERGED, TRIMMED_MERGED, TRIMMED_QC, MAPPED_BOWTIE2_PE_SAM, MAPPED_BOWTIE2_PE_BAM,
#     input: GENOMECOV, GENOMECOV_PLUS, GENOMECOV_MINUS
# #        MAPPED_PE_SORTED,  \
# #        GENOMECOV, \
# #        HTSEQ_COUNTS, \
# #        FEATUREOUNTS, \
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
# MAPPED_BOWTIE2_PE_SAM_OL=report_numbered_list(MAPPED_BOWTIE2_PE_SAM)
# MAPPED_BOWTIE2_PE_BAM_OL=report_numbered_list(MAPPED_BOWTIE2_PE_BAM)
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
#         {MAPPED_BOWTIE2_PE_SAM_OL}
#
#         Bam format (compressed)
#
#         {MAPPED_BOWTIE2_PE_BAM_OL}
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
