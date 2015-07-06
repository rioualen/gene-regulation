"""Workflow for the analysis of transcriptomic RNA-seq data in
Desulfovibrio desulfuricans.

Reference genome: Desulfovibrio desulfuricans ATCC_27774 uid59213

Collaboration with Alain Dolla
Data from : Jeff Cole

"""

#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
import subprocess
from snakemake.utils import R
configfile: "scripts/snakefiles/workflows/no-stress_analysis_config.json"
workdir: os.getcwd()
## A MODIFIER
#LIST_ALL_COUNTS = "data/1258-BRM/N1/N1_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/N2/N2_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/N4/N4_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/NN2/NN2_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/NN4/NN4_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/NN5/NN5_bowtie2_mm1_sorted_name_count.txt"



#================================================================#
#                         Includes                               #
#================================================================#

include: "fg-chip-seq/scripts/snakefiles/rules/util.py"
# include: "fg-chip-seq/scripts/snakefiles/rules/util.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/flowcharts.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/gunzip.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/rsync.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/fastqc.rules"
# include: "scripts/snakefiles/rules/fastqc.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/sickle_paired_ends.rules"
# include: "scripts/snakefiles/rules/trimming.rules"
# include: "fg-chip-seq/scripts/snakefiles/rules/bowtie_build.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/bowtie2_paired_ends.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/convert_sam_to_bam.rules"
# include: "scripts/snakefiles/rules/convert_sam_to_bam.rules"
include: "scripts/snakefiles/rules/sorted_bam.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/htseq.rules"
# include: "fg-chip-seq/scripts/snakefiles/rules/featurecounts.rules"
# include: "scripts/snakefiles/rules/count_table.rules"
include: "fg-chip-seq/scripts/snakefiles/rules/index_bam.rules"
# include: "scripts/snakefiles/rules/index_bam.rules"

#================================================================#
#     Global variables                                           #
#================================================================#

# Usage: snakemake -c "qsub {params.qsub}" -j 25 --allow-ambiguit
# workdir:"/home/desulfo-no" # Root directory for the project. Should be adapted for porting.
# workdir: "/home/desulfo-no/no-stress_project/data/1258-BRM" # Root directory for the project. Should be adapted for porting.
# workdir:"/home/no-stress/Documents/no-stress_project/data/1258-BRM" # Local Root directory for the project. Should be adapted for porting.
workdir: os.getcwd() # Local Root directory for the project. Should be adapted for porting.
SAMPLE_IDS=read_sample_ids(config["files"]["samples"], column=1, verbose=0)

#DATASETS_FWD = expand("{sample_id}" + config["suffix"]["reads_fwd"], sample_id=SAMPLE_IDS)
#DATASETS_REV = expand("{sample_id}" + config["suffix"]["reads_rev"], sample_id=SAMPLE_IDS)
#DATASETS = DATASETS_FWD + DATASETS_REV

DATASETS = "N1_1 N1_2 N2_1 N2_2 N4_1 N4_2 S1_1 S1_2 S4_1 S4_2 S5_1 S5_2 NN2_1 NN2_2 NN4_1 NN4_2 NN5_1 NN5_2 SN1_1 SN1_2 SN2_1 SN2_2 SN5_1 SN5_2".split() # list of files
DATA_DIRS = "N1 N1 N2 N2 N4 N4 S1 S1 S4 S4 S5 S5 NN2 NN2 NN4 NN4 NN5 NN5 SN1 SN1 SN2 SN2 SN5 SN5".split() # list all directories

RAWR_FILES_FWD, RAWR_DIRS_FWD, RAWR_BASENAMES_FWD = glob_multi_dir(SAMPLE_IDS, "*" + config["suffix"]["reads_fwd"] + ".fastq.gz", config["dir"]["data_root"], config["suffix"]["reads_fwd"] + ".fastq.gz")

RAWR_FILES, RAWR_DIRS, RAWR_BASENAMES = glob_multi_dir(SAMPLE_IDS, "*_?" + ".fastq.gz", config["dir"]["data_root"], ".fastq.gz")

# Variables for EdgeR (Differential expression) 
COND_1 = config["Diff_Exp"]["cond1"]
COND_2 = config["Diff_Exp"]["cond2"]

# COUNT_FILES = "data/1258-BRM/N1/N1_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/N2/N2_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/N4/N4_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/NN2/NN2_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/NN4/NN4_bowtie2_mm1_sorted_name_count.txt", "data/1258-BRM/NN5/NN5_bowtie2_mm1_sorted_name_count.txt"
# COUNT_FILES = expand("data/1258-BRM/{sample_id}/{sample_id}_bowtie2_mm1_sorted_name_count.txt", sample_id = SAMPLE_IDS)
COUNT_FILES = expand(config["dir"]["results"] + "{sample_id}/{sample_id}_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + "_count.txt", sample_id = SAMPLE_IDS)
COUNT_RESULTS_EDGER = expand(config["dir"]["results"] + "DEG/{cond_1}_vs_{cond_2}/{cond_1}_VS_{cond_2}_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + config["edgeR"]["software"] +".tab", zip, cond_1=COND_1, cond_2=COND_2)

RESULTS_DESEQ2 = expand(config["dir"]["results"] + "DEG/{cond_1}_vs_{cond_2}/{cond_1}_VS_{cond_2}_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + config["DESeq2"]["software"] + ".tab", zip, cond_1=COND_1, cond_2=COND_2)
LOGS_DESEQ2 = expand(config["dir"]["results"] + "DEG/{cond_1}_VS_{cond_2}_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + config["DESeq2"]["software"] + ".log", zip, cond_1=COND_1, cond_2=COND_2)

PARAMS_R = config["dir"]["results"] + "DEG/sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + "_params.R"
ALL_COUNTS = config["dir"]["results"] + "DEG/sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_pe_sorted_" + config["htseq"]["order"] + "_allcounts.tab"

include: "scripts/snakefiles/rules/HTseq_allcount_params.rules"
include: "scripts/snakefiles/rules/edgeR.rules"
include: "scripts/snakefiles/rules/DESeq2.rules"

#================================================================#
#                         Workflow                               #
#================================================================#

RAW_FASTQC=expand(config["dir"]["results"] + "{data_dir}/{dataset}_fastqc/", zip, dataset=DATASETS, data_dir=DATA_DIRS)
# RAW_FASTQC=expand(config["dir"]["data_root"] + "{data_dir}/{dataset}_fastqc" + config["fastqc"]["extension"] + "/", zip, dataset=DATASETS, data_dir=DATA_DIRS)




rule all:
    """Run workflow for each replica of each experience"""
    input: RAW_FASTQC, \
        expand(config["dir"]["results"] + "{data_dir}/{dataset}_sickle_pe_q" + config["sickle"]["threshold"] + "_fastqc/", zip, dataset=DATASETS, data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_trimmed_thr" + THRESHOLD + ".fastq.gz", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + config["bowtie2"]["max_mismatches"] + ".sam", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + config["bowtie2"]["max_mismatches"] + ".bam", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + config["bowtie2"]["max_mismatches"] + "_sorted_" + config["htseq"]["order"] + ".bam", data_dir=DATA_DIRS), \
        # expand(config["dir"]["results"] + "{data_dir}/{data_dir}_sickle_pe_q" + config["sickle"]["threshold"] + "_bowtie2_sorted_pos.bam.bai", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + config["bowtie2"]["max_mismatches"] + "_sorted_" + config["htseq"]["order"] + "_count.txt", data_dir=DATA_DIRS)
        PARAMS_R, \
        ALL_COUNTS, \
        COUNT_RESULTS_EDGER, \
        RESULTS_DESEQ2


    shell: "echo 'job done'"



# ruleorder: bowtie2_paired_end > sorted_bam
# ruleorder: sorted_bam > index_bam

ruleorder: rsync > sickle_paired_ends




################################################################
## Build the report (including DAG and rulegraph flowcharts).
from snakemake.utils import report

## Bulleted list of samples for the report
SAMPLE_IDS_OL=report_numbered_list(SAMPLE_IDS)
RAW_FASTQC_OL=report_numbered_list(RAW_FASTQC)

rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  dag_pdf=config["reports"]["dag"] + ".pdf", \
            dag_png=config["reports"]["dag"] + ".png", \
            rulegraph_pdf=config["reports"]["rulegraph"] + ".pdf", \
            rulegraph_png=config["reports"]["rulegraph"] + ".png"
    output: html=config["dir"]["reports"] + "/report.html"
    run:
        report("""
        ============================
        RNA-seq analysis - Jeff Cole
        ============================
        
        :Project:              Jeff Cole
        :Collaboration:        Alain Dolla
        :Analysis workflow:    Justine Long, Jeanne Chèneby & Jacques van Helden
        
        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets`_
             - `Sample identifiers`_

        -----------------------------------------------------

        Flowcharts
        ==========

        - Sample treatment: dag_pdf_
        - Workflow: rulegraph_pdf_

        .. image:: rulegraph.png

        -----------------------------------------------------

        Datasets
        ========
        
        Sample identifiers
        ------------------

        {SAMPLE_IDS_OL} 

        Result files
        ============

        Quality control (raw reads)
        ---------------------------

        {RAW_FASTQC_OL}

        Raw reads (forward)
        -------------------
        {RAWR_FILES_FWD}

        Raw read directories (forward)
        ------------------------------
        {RAWR_DIRS_FWD}

        Raw read basenames (forward)
        ----------------------------
        {RAWR_BASENAMES_FWD}

        Raw reads 
        -------------------
        {RAWR_FILES}

        Raw read directories
        ------------------------------
        {RAWR_DIRS}

        Raw read basenames
        ----------------------------
        {RAWR_BASENAMES}

        """, output.html, metadata="Jacques van Helden (Jacques.van-Helden@univ-amu.fr)", **input)
