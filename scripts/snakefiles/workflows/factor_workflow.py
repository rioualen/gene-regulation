"""Generic work flow for the analysis of ChIP-seq data for the binding
of transcription factors.


This workflow performs the following treatments: 

 - read mapping
 - peak-calling with alternate peak-calling programs
 - motif discovery

The details are specified in a yaml-formatted configuration file.

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/factor_workflow.py \
        --configfile path/to/specific/config_file.yml \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/factor_workflow.py \
        --configfile path/to/specific/config_file.yml \
        --force flowcharts

Reference genome:	-
Sequencing type: 	single end

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
configfile: "examples/GSE20870/GSE20870.yml"
#workdir: config["dir"]["base"]
#verbosity = int(config["verbosity"])

#================================================================#
#                    Check mandatory parameters
#================================================================#

# Define verbosity
if not ("verbosity" in config.keys()):
    config["verbosity"] = 0
verbosity = int(config["verbosity"])

#================================================================#
#                         Includes                               #
#================================================================#
if not ("dir" in config.keys()) & ("fg_lib" in config["dir"].keys()) :
    sys.exit("The parameter config['dir']['fg_lib'] should be specified in the config file.")
FG_LIB = os.path.abspath(config["dir"]["fg_lib"])
RULES = os.path.join(FG_LIB, "scripts/snakefiles/rules")
PYTHON = os.path.join(FG_LIB, "scripts/python_lib")

include: os.path.join(PYTHON, "util.py")

include: os.path.join(RULES, "bam_by_name.rules")
include: os.path.join(RULES, "bam_by_pos.rules")
include: os.path.join(RULES, "bam_to_bed.rules")
include: os.path.join(RULES, "bowtie2_index.rules")
include: os.path.join(RULES, "bowtie2_se.rules")
include: os.path.join(RULES, "bPeaks.rules")
include: os.path.join(RULES, "bwa_index.rules")
include: os.path.join(RULES, "bwa_se.rules")
include: os.path.join(RULES, "count_reads.rules")
#include: os.path.join(RULES, "download_from_GEO.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "flowcharts.rules")
include: os.path.join(RULES, "genome_download.rules")
include: os.path.join(RULES, "getfasta.rules")
include: os.path.join(RULES, "get_chrom_sizes.rules")
#include: os.path.join(RULES, "genome_coverage.rules")
include: os.path.join(RULES, "homer.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "macs14.rules")
include: os.path.join(RULES, "merge_lanes.rules")               ## Merge lanes by sample, based on a tab-delimited file indicainclude: os.path.join(RULES, "peak_motifs.rules")
include: os.path.join(RULES, "peak_motifs.rules")
include: os.path.join(RULES, "sra_to_fastq.rules")
include: os.path.join(RULES, "swembl.rules")
include: os.path.join(RULES, "sam_to_bam.rules")

ruleorder: bam_by_pos > sam_to_bam
ruleorder: bam_by_name > sam_to_bam

#================================================================#
#                      Data & wildcards                          #
#================================================================#

# Raw data
READS = config["dir"]["reads_source"]

# Samples
SAMPLES = read_table(config["files"]["samples"])
SAMPLE_IDS = SAMPLES.iloc[:,0] ## First column MUST contain the sample ID
SRR_IDS = SAMPLES['SRR']

## Design
DESIGN = read_table(config["files"]["design"])
TREATMENT = DESIGN['treatment']
CONTROL = DESIGN['control']

## Genome
#GENOME = config["genome"]["version"]
#GENOME_DIR = config["genome"]["dir"]
#CHROM_SIZES = expand(GENOME_DIR + "{genome}.genome", genome=GENOME)
#config["genome"]["chromsize"] = CHROM_SIZES

## Data & results dir

if not (("dir" in config.keys()) and ("reads_source" in config["dir"].keys())):
    sys.exit("The parameter config['dir']['reads_source'] should be specified in the config file.")

READS = config["dir"]["reads_source"]
if not os.path.exists(READS):
    os.makedirs(READS)

if not ("results" in config["dir"].keys()):
    sys.exit("The parameter config['dir']['results'] should be specified in the config file.")

RESULTS_DIR = config["dir"]["results"]
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

if not ("samples" in config["dir"].keys()):
    SAMPLE_DIR = config["dir"]["results"]
else:
    SAMPLE_DIR = config["dir"]["samples"]
if not os.path.exists(SAMPLE_DIR):
    os.makedirs(SAMPLE_DIR)

if not ("peaks" in config["dir"].keys()):
    PEAKS_DIR = config["dir"]["results"]
else:
    PEAKS_DIR = config["dir"]["peaks"]
if not os.path.exists(PEAKS_DIR):
    os.makedirs(PEAKS_DIR)


#================================================================#
#                         Workflow                               #
#================================================================#

## Data import & merging.


#DOWNLOAD = expand(READS + "{samples}/{srr}.sra", zip, samples=SAMPLE_IDS, srr=SRR_IDS)
IMPORT = expand(SAMPLE_DIR + "{samples}/{samples}.fastq", zip, samples=SAMPLE_IDS, srr=SRR_IDS)


## Verbosity
#if (verbosity >= 3):
#    print("\tDOWNLOAD:\t" + ";".join(DOWNLOAD))
#    print("\tIMPORT:\t" + ";".join(IMPORT))

### Graphics & reports
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")
#REPORT = expand(RESULTS_DIR + "report.html")

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(SAMPLE_DIR + "{samples}/{samples}_fastqc/{samples}_fastqc.html", samples=SAMPLE_IDS)
RAW_READNB = expand(SAMPLE_DIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)

#----------------------------------------------------------------#
# Alignment
#----------------------------------------------------------------#


ALIGNER=["bowtie2", "bwa"]
ALIGNMENT=expand(SAMPLE_DIR + "{samples}/{samples}_{aligner}", samples=SAMPLE_IDS, aligner=ALIGNER)


INDEX = expand(config["dir"]["genomes"] + config["genome"]["version"] + "/{aligner}/" + config["genome"]["version"] + ".fa", aligner=ALIGNER)

#BWA_INDEX = expand(config["dir"]["genomes"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
#BOWTIE2_INDEX = expand(config["dir"]["genomes"] + "{genome}/Bowtie2Index/{genome}.fa.1.bt2", genome=GENOME)

MAPPING = expand("{alignment}.sam", alignment=ALIGNMENT)


# Sorted and converted reads (bam, bed)
SORTED_MAPPED_READS_BWA = expand(SAMPLE_DIR + "{alignment}_sorted_pos.bam", alignment=ALIGNMENT)
#BAM_READNB = expand(RESULTS_DIR + "{alignment}_sorted_pos_bam_readnb.txt", alignment=ALIGNMENT)
SORTED_READS_BED = expand(SAMPLE_DIR + "{alignment}_sorted_pos.bed", alignment=ALIGNMENT)
#BED_FEAT_COUNT = expand(RESULTS_DIR + "{alignment}_sorted_pos_bed_nb.txt", alignment=ALIGNMENT)

#TDF = expand(RESULTS_DIR + "{alignment}_sorted_pos.tdf", alignment=ALIGNMENT)

# ----------------------------------------------------------------
# Peak-calling
# ----------------------------------------------------------------



PEAKCALLER=[
    "homer-fdr" + config["homer"]["fdr"] + "_peaks", 
    "macs2-qval" + config["macs2"]["qval"], 
    "swembl-R" + config["swembl"]["R"],
    "macs14-pval" + config["macs14"]["pval"],
#    "spp-fdr" + config["spp"]["fdr"],
    "bPeaks_allGenome"
]

PEAKCALLING=expand(expand(PEAKS_DIR + "{treat}_vs_{control}/{{peakcaller}}/{treat}_vs_{control}_{{aligner}}_{{peakcaller}}", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER)

PEAKS = expand("{peakcalling}.bed", peakcalling=PEAKCALLING)

# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------


MOTIFS=expand(expand("{treat}_vs_{control}/{{peakcaller}}/peak-motifs/{treat}_vs_{control}_{{aligner}}_{{peakcaller}}_peak-motifs_synthesis", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER)


GET_FASTA = expand(PEAKS_DIR + "{peakcalling}.fasta", peakcalling=PEAKCALLING)
PEAK_MOTIFS = expand(PEAKS_DIR + "{motifs}.html", motifs=MOTIFS)

#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses.
	"""
	input: IMPORT, INDEX, MAPPING#, PEAKS, PEAK_MOTIFS, GRAPHICS#, CHROM_SIZES, PEAKS, TDFRAW_QC, 
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"
