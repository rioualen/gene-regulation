"""This workflow was designed in order to build and test a pipeline
based on histone marks ChIP-seq analysis.  Data was downloaded from
the GEO platform. It is composed of 2 inputs and 2 ChIP targetting
H3K27me3 in C. elegans.

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/Paeruginosa.py \
        --force flowcharts

Organism: 		Saccharomyces cerevisiae
Reference genome:	sacCer3
Sequencing type: 	single end
Data source: 		Gene Expression Omnibus

Author: 		Claire Rioualen, Jacques van Helden
Contact: 		claire.rioualen@inserm.fr
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

## Config
configfile: "scripts/snakefiles/workflows/Scerevisiae_new.json"
workdir: config["dir"]["base"]
verbosity = int(config["verbosity"])

#================================================================#
#                         Includes                               #
#================================================================#

FG_LIB = os.path.abspath(config["dir"]["fg_lib"])
RULES = os.path.join(FG_LIB, "scripts/snakefiles/rules")
PYTHON = os.path.join(FG_LIB, "scripts/snakefiles/python_lib")

include: os.path.join(PYTHON, "util.py")
include: os.path.join(RULES, "util.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "bwa_index.rules")
include: os.path.join(RULES, "bwa_se.rules")
#include: os.path.join(RULES, "bowtie2_index.rules")
#include: os.path.join(RULES, "bowtie2_se.rules")
include: os.path.join(RULES, "convert_bam_to_bed.rules")
include: os.path.join(RULES, "count_oligo.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "flowcharts.rules")
include: os.path.join(RULES, "getfasta.rules")
include: os.path.join(RULES, "homer.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "peak_length.rules")
include: os.path.join(RULES, "peak_motifs.rules")
include: os.path.join(RULES, "purge_sequence.rules")
include: os.path.join(RULES, "sickle_se.rules")
include: os.path.join(RULES, "spp.rules")
include: os.path.join(RULES, "sra_to_fastq.rules")
include: os.path.join(RULES, "swembl.rules")

#================================================================#
#                      Data & config                             #
#================================================================#

# Raw data
READS = config["dir"]["reads_source"]

# Samples
SAMPLES = read_table(config["files"]["samples"], verbosity=verbosity)
SAMPLE_IDS = SAMPLES.iloc[:,0] ## First column MUST contain the sample ID

## Design
DESIGN = read_table(config["files"]["design"], verbosity=verbosity)
TREATMENT = DESIGN.iloc[:,0]
CONTROL = DESIGN.iloc[:,1]

## Ref genome
GENOME = config["genome"]["version"]

## Results dir
RESULTS_DIR = config["dir"]["results"]
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

#================================================================#
#                         Workflow                               #
#================================================================#

## Data import & merging.

IMPORT = expand(RESULTS_DIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS) 

## Graphics & reports
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")
REPORT = expand(RESULTS_DIR + "report.html")

## Suffixes (beta) (! implement several values for each param)

TRIMMER="sickle-se-q" + config["sickle"]["threshold"]
TRIMMING=expand("{samples}/{samples}_{trimmer}", samples=SAMPLE_IDS, trimmer=TRIMMER)

ALIGNER="bwa"
ALIGNMENT=expand("{samples}/{samples}_{trimmer}_{aligner}", samples=SAMPLE_IDS, trimmer=TRIMMER, aligner=ALIGNER)

PEAKCALLER="homer_peaks swembl-R" + config["swembl"]["R"] # macs2_peaks spp-fdr" + config["spp"]["fdr"] + " 
PEAKCALLER=PEAKCALLER.split()
PEAKCALLING=expand(expand("{treat}_vs_{control}/{{peakcaller}}/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, trimmer=TRIMMER, aligner=ALIGNER)

MOTIFS=expand(expand("{treat}_vs_{control}/{{peakcaller}}/peak-motifs/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}_purged", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, trimmer=TRIMMER, aligner=ALIGNER)

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(RESULTS_DIR + "{samples}/{samples}_fastqc/", samples=SAMPLE_IDS)
RAW_READNB = expand(RESULTS_DIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)

TRIM = expand(RESULTS_DIR + "{trimming}.fastq", trimming=TRIMMING)
TRIM_QC = expand(RESULTS_DIR + "{samples}/{samples}_{trimmer}_fastqc/", samples=SAMPLE_IDS, trimmer=TRIMMER)

#----------------------------------------------------------------#
# Alignment
#----------------------------------------------------------------#

## to avoid duplicates, fasta sequence should be moved to {genome} directly...
BWA_INDEX = expand(config["dir"]["genome"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
#BOWTIE2_INDEX = expand(config["dir"]["genome"] + "{genome}/Bowtie2Index/{genome}.fa.1.bt2", genome=GENOME)

MAPPING = expand(RESULTS_DIR + "{alignment}.sam", alignment=ALIGNMENT)

# Sorted and converted reads (bam, bed)
SORTED_MAPPED_READS_BWA = expand(RESULTS_DIR + "{alignment}_sorted_pos.bam", alignment=ALIGNMENT)
BAM_READNB = expand(RESULTS_DIR + "{alignment}_sorted_pos_bam_readnb.txt", alignment=ALIGNMENT)
SORTED_READS_BED = expand(RESULTS_DIR + "{alignment}_sorted_pos.bed", alignment=ALIGNMENT)
BED_FEAT_COUNT = expand(RESULTS_DIR + "{alignment}_sorted_pos_bed_nb.txt", alignment=ALIGNMENT)

# ----------------------------------------------------------------
# Peak-calling
# ----------------------------------------------------------------

#BEDS = expand(RESULTS_DIR + "{alignment}.bed", alignment=ALIGNMENT)
PEAKS = expand(RESULTS_DIR + "{peakcalling}.bed", peakcalling=PEAKCALLING)

# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------

GET_FASTA = expand(RESULTS_DIR + "{peakcalling}.fasta", peakcalling=PEAKCALLING)
PURGE_PEAKS = expand(RESULTS_DIR + "{peakcalling}_purged.fasta", peakcalling=PEAKCALLING)
PEAKS_LENGTH = expand(RESULTS_DIR + "{peakcalling}_purged_length.png", peakcalling=PEAKCALLING)
PEAK_MOTIFS = expand(RESULTS_DIR + "{motifs}_peak-motifs_synthesis.html", motifs=MOTIFS)

## Oligo analysis # ! missing f* input exception
OLIGO = config['oligo_analysis']['count_oligo'].split()
OLIGO_ANALYSIS = expand(RESULTS_DIR + "{peakcalling}_purged_oligo{oligo}.txt", peakcalling=PEAKCALLING, oligo=OLIGO)

#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses
	"""
	input: GRAPHICS, RAW_QC, TRIM_QC, PEAK_MOTIFS
	#BED_FEAT_COUNT, PURGE_PEAKS, PEAKS_LENGTH
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

#================================================================#
#                          Report                                #
#================================================================#

NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

#rule report:
#    """
#    Generate a report with the list of datasets + summary of the results.
#    """
# see Scerevisiae report

#----------------------------------------------------------------#
# Build the report (including DAG and rulegraph flowcharts).
from snakemake.utils import report

# Bulleted list of samples for the report
SAMPLE_IDS_OL=report_numbered_list(SAMPLE_IDS)
RAW_READS_OL=report_numbered_list(IMPORT)
TRIMMED_READS_OL=report_numbered_list(TRIMMED_READS_SICKLE)
RAW_QC_OL=report_numbered_list(RAW_QC)
TRIMMED_QC_OL=report_numbered_list(TRIMMED_QC)

MAPPED_SAM_OL=report_numbered_list(MAPPED_READS_BWA)
MAPPED_BAM_SORTED_OL=report_numbered_list(SORTED_MAPPED_READS_BWA)
MAPPED_BED_SORTED_OL=report_numbered_list(SORTED_READS_BED)

RAW_QC_OL=report_numbered_list(RAW_QC)
TRIMMED_QC_OL=report_numbered_list(TRIMMED_QC)

PEAKFILES_MACS2_OL=report_numbered_list(PEAKS_MACS2)

#	input: GRAPHICS, IMPORT, TRIMMED_READS_SICKLE, TRIMMED_QC, RAW_QC, MAPPED_READS_BWA, RAW_READNB, BAM_READNB, BED_READNB, PEAKS_MACS2, FETCH_MACS2_PEAKS, PURGE_MACS2_PEAKS #redundant for flowcharts

NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

rule report:
    """
    Generate a report with the list of datasets + summary of the results.
    """
    input:  dag=config["dir"]["reports"] + "dag.pdf", \
            dag_png=config["dir"]["reports"] + "dag.png", \
            rulegraph=config["dir"]["reports"] + "rule.pdf", \
            rulegraph_png=config["dir"]["reports"] + "rule.png"
    output: html=config["dir"]["reports"] + "report.html"
    run:
        report("""
        ===========================================
        ChIP-seq analysis - S.cerevisiae Pho4 study
        ===========================================
        
        :Date:                 {NOW}
        :Project:              S cerevisiae
        :Analysis workflow:    Claire Rioualen
        
        Contents
        ========
        
        - `Flowcharts`_
        - `Datasets`_
             - `Samples`_
             - `Raw reads`_
             - `Trimmed`_
             - `Mapped`_
             - `Peaks`_
             - `QC reports`_

        -----------------------------------------------------

        Flowcharts
        ==========

        - Sample treatment: dag_
        - Workflow: rulegraph_

        .. image:: rulegraph.png

        -----------------------------------------------------

        Datasets
        ========
        
        Samples
        -------

        {SAMPLE_IDS_OL} 

        Raw reads 
        ---------

        {RAW_READS_OL}

        Trimmed
        -------

        {TRIMMED_READS_OL}

        Mapped
        ------

        Sam format (uncompressed)

        {MAPPED_SAM_OL}

        Bam format (sorted by positions)

        {MAPPED_BAM_SORTED_OL}

        Bed format (sorted by positions)

        {MAPPED_BED_SORTED_OL}

        Peaks
        -----

        Macs2 peaks

        {PEAKFILES_MACS2_OL}

        QC reports
        ----------

        {RAW_QC_OL}

        {TRIMMED_QC_OL}

        -----------------------------------------------------

        """, output.html, metadata="Claire Rioualen (claire.rioualen@inserm.fr)", **input)




