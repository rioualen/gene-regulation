"""This workflow was designed in order to build and test a pipeline
based on histone marks ChIP-seq analysis.  Data was downloaded from
the GEO platform. It is composed of 2 inputs and 2 ChIP targetting
H3K27me3 in C. elegans.

The workflows should include: 
	- normalization & QC: fastq merging, trimming, quality control...
	- formatting rules: bed, bam, sam, narrowPeak, fasta...
	- alignment: Bowtie, BWA
	- peak-calling: SWEMBL, MACS, SPP, HOMER
	- downstream analyses: IDR, sequence purge, peak length

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        --force flowcharts

Organism: 		Saccharomyces cerevisiae
Reference genome:	
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

#================================================================#
#                      Data & directories                        #
#================================================================#

configfile: "scripts/snakefiles/workflows/Scerevisiae-Pho4_VM.json"
workdir: config["dir"]["base"]

# Beware: verbosity messages are incompatible with the flowcharts
verbosity = int(config["verbosity"])

# Raw data
READS = config["dir"]["reads_source"]
RESULTS_DIR = config["dir"]["results"]	
if verbosity >= 1:
    print("\nRESULTS_DIR\t" + RESULTS_DIR)

GENOME = config["genome"]["genome_version"]
if verbosity >= 1:
    print("\nGENOME\t" + GENOME)

#================================================================#
#                         Includes                               #
#================================================================#

if verbosity >= 2:
    print("\nImporting rules")

#if "base" in config["dir"]:
#    workflow.basedir = config["dir"]["base"]
#    print("workflow.basedir\t" + workflow.basedir)
#    print('srcdir("broche_analysis.py")' + "\t" + srcdir("broche_analysis.py"))
#    print('srcdir("../python_lib/util.py")' + "\t" + srcdir("../python_lib/util.py"))


FG_LIB = os.path.abspath(config["dir"]["fg_lib"])
RULES = os.path.join(FG_LIB, "scripts/snakefiles/rules")
PYTHON = os.path.join(FG_LIB, "scripts/snakefiles/python_lib")
#PYTHON = config["dir"]["python_lib"]
#PYTHON = "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/python_lib/"

#print ("PYTHON\t" + PYTHON)
#print('config["dir"]["base"]' + "\t" + config["dir"]["base"])
#print("CWD\t" + os.getcwd())

include: os.path.join(PYTHON, "util.py")
include: os.path.join(RULES, "util.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "bwa_index.rules")
include: os.path.join(RULES, "bwa_se.rules")
include: os.path.join(RULES, "bed_to_fasta.rules")
include: os.path.join(RULES, "convert_bam_to_bed.rules")
include: os.path.join(RULES, "count_oligo.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "flowcharts.rules")
include: os.path.join(RULES, "homer.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "peak_length.rules")
include: os.path.join(RULES, "peak_motifs.rules")
include: os.path.join(RULES, "purge_sequence.rules")
include: os.path.join(RULES, "sickle_se.rules")
#include: os.path.join(RULES, "sorted_bam.rules")
include: os.path.join(RULES, "spp.rules")
include: os.path.join(RULES, "sra_to_fastq.rules")
include: os.path.join(RULES, "swembl.rules")


#================================================================#
#                Read sample IDs and design                      #
#================================================================#


## Read sample descriptions
if verbosity >= 1:
    print("\nReading sample descriptions\t" + config["files"]["samples"])

# Read the sample description file
SAMPLES = read_table(config["files"]["samples"], verbosity=verbosity)
SAMPLE_IDS = SAMPLES.iloc[:,0] ## First column MUST contain the sample ID
SAMPLE_CONDITIONS = SAMPLES.iloc[:,1] ## Second column MUST contain condition for each sample
SAMPLE_EXT_DB = SAMPLES.iloc[:,2] ## External database from which the raw reads will be downloaded
SAMPLE_EXT_ID = SAMPLES.iloc[:,3] ## Sample ID in the external database
SAMPLE_DESCR = SAMPLES.iloc[:,7] ## First column MUST contain the sample ID

if verbosity >= 2:
    print("\nSAMPLE_IDS\t" + ";".join(SAMPLE_IDS))

## Read the list of sample IDs from the sample description file
if verbosity >= 1:
    print("\nReading experimental design\t" + config["files"]["design"])

DESIGN = read_table(config["files"]["design"], verbosity=verbosity)
TREATMENT = DESIGN.iloc[:,0]
CONTROL = DESIGN.iloc[:,1]

if verbosity >= 2:
    print("\nTREATMENT\t" + ";".join(TREATMENT))
    print("\nCONTROL\t" + ";".join(CONTROL))

#================================================================#
#                         Workflow                               #
#================================================================#

# ## Data import & merging.

## TODO
# Moved to sra_to_fastq.rules
# Should develop general import rule?
IMPORT = expand(RESULTS_DIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS) 
	

# Graphics & reports
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")
REPORT = expand(RESULTS_DIR + "report.html")

#----------------------------------------------------------------#
# Trimming of the raw reads
#----------------------------------------------------------------#

TRIMMED_READS_SICKLE = expand(RESULTS_DIR + "{samples}/{samples}" + config["sickle"]["suffix"] + ".fastq", samples=SAMPLE_IDS)
if verbosity >= 3:
    print("\nTRIMMED_READS_SICKLE\n\t" + "\n\t".join(TRIMMED_READS_SICKLE))

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(RESULTS_DIR + "{samples}/{samples}_fastqc/", samples=SAMPLE_IDS)
if verbosity >= 3:
    print("\nRAW_QC\n\t" + "\n\t".join(RAW_QC))
RAW_READNB = expand(RESULTS_DIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)

TRIMMED_QC = expand(RESULTS_DIR + "{samples}/{samples}" + config["sickle"]["suffix"] + "_fastqc/", samples=SAMPLE_IDS)
if verbosity >= 3:
    print("\nTRIMMED_QC\n\t" + "\n\t".join(TRIMMED_QC))

#----------------------------------------------------------------#
# Mapping
#----------------------------------------------------------------#

BWA_INDEX = expand(config["dir"]["genome"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
MAPPED_READS_BWA = expand(RESULTS_DIR + "{samples}/{samples}" + config["sickle"]["suffix"] + "_{aligner}.sam", samples=SAMPLE_IDS, aligner=config["bwa"]["suffix"])
if verbosity >= 3:
    print("\nMAPPED_READS_BWA\n\t" + "\n\t".join(MAPPED_READS_BWA))

# Sorted and converted reads (bam, bed)
SORTED_MAPPED_READS_BWA = expand(RESULTS_DIR + "{sample}/{sample}" + config["sickle"]["suffix"] + "_{aligner}_sorted_pos.bam", sample=SAMPLE_IDS, aligner=config["bwa"]["suffix"])
BAM_READNB = expand(RESULTS_DIR + "{samples}/{samples}" + config["sickle"]["suffix"] + "_{aligner}_sorted_pos_bam_readnb.txt", samples=SAMPLE_IDS, aligner=config["bwa"]["suffix"])
SORTED_READS_BED = expand(RESULTS_DIR + "{sample}/{sample}" + config["sickle"]["suffix"] + "_{aligner}_sorted_pos.bed", sample=SAMPLE_IDS, aligner=config["bwa"]["suffix"])
BED_READNB = expand(RESULTS_DIR + "{samples}/{samples}" + config["sickle"]["suffix"] + "_{aligner}_sorted_pos_bed_nb.txt", samples=SAMPLE_IDS, aligner=config["bwa"]["suffix"])
if verbosity >= 3: 
    print("\nSORTED_READS_BED\n\t" + "\n\t".join(SORTED_READS_BED))

# ----------------------------------------------------------------
# Peak-calling
# ----------------------------------------------------------------

## Peak-calling with MACS2
PEAKS_MACS2 = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks.narrowPeak",
               zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
if verbosity >= 3: 
    print("\nPEAKS_MACS2\n\t" + "\n\t".join(PEAKS_MACS2))

## Peak-calling with SWEMBL
PEAKS_SWEMBL = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/swembl/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_swembl-R" + config["swembl"]["R"] + ".bed", 
               zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
if verbosity >= 3: 
    print("\nPEAKS_SWEMBL\n\t" + "\n\t".join(PEAKS_SWEMBL))

## Peak-calling with SPP
PEAKS_SPP = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/spp/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_spp.narrowPeak", 
               zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
if verbosity >= 3: 
    print("\nPEAKS_SPP\n\t" + "\n\t".join(PEAKS_SPP))

PEAKS_HOMER = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks.bed", 
               zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
if verbosity >= 3: 
    print("\nPEAKS_SPP\n\t" + "\n\t".join(PEAKS_SPP))

## Peaks returned by the different peak callers
PEAKS = PEAKS_MACS2 + PEAKS_HOMER + PEAKS_SWEMBL

# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------

# File conversion / fetching fasta
FETCH_MACS2_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
FETCH_SWEMBL_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/swembl/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_swembl-R0.01.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
FETCH_SPP_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/spp/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_spp.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
FETCH_HOMER_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])

FETCH_PEAKS = FETCH_MACS2_PEAKS + FETCH_HOMER_PEAKS

# Sequence purge
PURGE_MACS2_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
PURGE_SWEMBL_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/swembl/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_swembl-R0.01_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
PURGE_SPP_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/spp/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_spp_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
PURGE_HOMER_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])

PURGE_PEAKS = PURGE_MACS2_PEAKS + PURGE_HOMER_PEAKS

## Working til here ##

# Oligo analysis ## BUG oligo.sam...
OLIGO = config['oligo_analysis']['count_oligo'].split()
MACS2_OLIGO = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged_oligo{{oligo}}.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"], oligo=OLIGO)
HOMER_OLIGO = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged_oligo{{oligo}}.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"], oligo=OLIGO)

# Peaks length ## TODO when purge works -> vmatch issue
MACS2_PEAKS_LENGTH = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged_length.png", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
HOMER_PEAKS_LENGTH = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged_length.png", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])

# Peak-motifs
MACS2_MOTIFS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged.html", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])
HOMER_MOTIFS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged.html", zip, treat=TREATMENT, ctrl=CONTROL), trimming=config["sickle"]["suffix"], aligner=config["bwa"]["suffix"])


ruleorder: macs2 > bam_to_bed > sam2bam
ruleorder: count_reads_bam > sam2bam
ruleorder: bed_to_fasta > purge_sequence

# tests ruleorder
ruleorder: swembl > bam_to_bed > sort_bam_by_pos > sam2bam


rule all: 
    """Run all the required analyses.

    """
#	input: GRAPHICS, IMPORT, TRIMMED_READS_SICKLE, TRIMMED_QC, RAW_QC, MAPPED_READS_BWA, RAW_READNB, BAM_READNB, BED_READNB, PEAKS_MACS2, FETCH_MACS2_PEAKS, PURGE_MACS2_PEAKS #redundant for flowcharts
	input: GRAPHICS, REPORT, TRIMMED_QC, RAW_QC, RAW_READNB, BAM_READNB, BED_READNB, PEAKS, PURGE_PEAKS
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

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




