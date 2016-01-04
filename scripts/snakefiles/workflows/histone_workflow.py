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

Organism: 		C elegans
Reference genome:	ce10 (WBcel235?)
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
#configfile: "examples/Celegans-H3K27me3/Celegans-H3K27me3.yml"
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
include: os.path.join(RULES, "macs14.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "peak_length.rules")
include: os.path.join(RULES, "peak_motifs.rules")
#include: os.path.join(RULES, "purge_sequence.rules")
include: os.path.join(RULES, "sickle_se.rules")
include: os.path.join(RULES, "sicer.rules")
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

# /!\ C.elegans data: chromosomes not properly named. Becomes an issue with fetch-sequences or IGV visualization for instance. 

#find . -type f -print0 | xargs -0 sed -i 's/NC_003279.8/chrI/g;s/NC_003280.10/chrII/g;s/NC_003281.10/chrIII/g;s/NC_003282.8/chrIV/g;s/NC_003283.11/chrV/g;s/NC_003284.9/chrX/g'

## Graphics & reports
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")
REPORT = expand(RESULTS_DIR + "report.html")


## Programs

TRIMMER="sickle-se-q" + config["sickle"]["threshold"]
TRIMMING=expand("{samples}/{samples}_{trimmer}", samples=SAMPLE_IDS, trimmer=TRIMMER)

ALIGNER="bwa"
ALIGNMENT=expand("{samples}/{samples}_{trimmer}_{aligner}", samples=SAMPLE_IDS, trimmer=TRIMMER, aligner=ALIGNER)

PEAKCALLER=[
    "homer_peaks", 
    "macs2-qval" + config["macs2"]["qval"] + "_peaks",
    "macs14-pval" + config["macs14"]["pval"] + "_peaks",#
    "swembl-R" + config["swembl"]["R"],
#    "sicer-fdr" + config["sicer"]["fdr"]
#    "spp-fdr" + config["spp"]["fdr"] todo option broadpeaks
]
PEAKCALLING=expand(expand("{treat}_vs_{control}/{{peakcaller}}/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER, trimmer=TRIMMER)

MOTIFS=expand(expand("{treat}_vs_{control}/{{peakcaller}}/peak-motifs/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, trimmer=TRIMMER, aligner=ALIGNER)

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(RESULTS_DIR + "{samples}/{samples}_fastqc/", samples=SAMPLE_IDS)
RAW_READNB = expand(RESULTS_DIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)

TRIM = expand(RESULTS_DIR + "{trimming}.fastq", trimming=TRIMMING)
TRIM_QC = expand(RESULTS_DIR + "{samples}/{samples}_{trimmer}_fastqc/", samples=SAMPLE_IDS, trimmer=TRIMMER)

QC = RAW_QC + TRIM_QC

#----------------------------------------------------------------#
# Alignment
#----------------------------------------------------------------#

BWA_INDEX = expand(config["genome"]["dir"] + "BWAIndex/{genome}.fa.bwt", genome=GENOME)
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
#PURGE_PEAKS = expand(RESULTS_DIR + "{peakcalling}_purged.fasta", peakcalling=PEAKCALLING)
#PEAKS_LENGTH = expand(RESULTS_DIR + "{peakcalling}_purged_length.png", peakcalling=PEAKCALLING)
PEAK_MOTIFS = expand(RESULTS_DIR + "{motifs}_peak-motifs_synthesis.html", motifs=MOTIFS)

## Oligo analysis # ! missing f* input exception
#OLIGO = config['oligo_analysis']['count_oligo'].split()
#OLIGO_ANA/data/results/Celegans-H3K4me3/GSM1208361_vs_GSM1217457/macs2-qval0.05_peaks/GSM1208361_vs_GSM1217457_sickle-se-q20_bwa_macs2-qval0.05.logLYSIS = expand(RESULTS_DIR + "{peakcalling}_purged_oligo{oligo}.txt", peakcalling=PEAKCALLING, oligo=OLIGO)

#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses
	"""
	input: GRAPHICS, QC, PEAKS#, PEAK_MOTIFS#, TRIM_QC IMPORT RAW_QC, TRIM_QC,   
	#BED_FEAT_COUNT, PURGE_PEAKS, PEAKS_LENGTH
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"


#ruleorder: chr_names > sam_to_bam

#rule chr_names:
#    input: "{file}.sam"
#    output: "{file}.bam"
#    shell:"awk '{{gsub(\"NC_003279.8\",\"chrI\"); gsub(\"NC_003280.10\",\"chrII\"); gsub(\"NC_003281.10\",\"chrIII\"); gsub(\"NC_003282.8\",\"chrIV\"); gsub(\"NC_003283.11\",\"chrV\"); gsub(\"NC_003284.9\",\"chrX\"); print}}' {input} > {input}.converted; \
#            samtools view -b -S {input}.converted > {output}"


## /!\ C.elegans chromosomes not properly named, this should be acted upon at the beginning of the workflow /!\
# tmp
#rule bed_chr_id_conversion:
#	input: "{file}.aligned.bed"
#	output:	"{file}.converted.bed"
#	shell:"awk '{{gsub(\"NC_003279.8\",\"chrI\"); gsub(\"NC_003280.10\",\"chrII\"); gsub(\"NC_003281.10\",\"chrIII\"); gsub(\"NC_003282.8\",\"chrIV\"); gsub(\"NC_003283.11\",\"chrV\"); gsub(\"NC_003284.9\",\"chrX\"); print}}' {input} > {output[0]}"

#find . -type f -print0 | xargs -0 sed -i 's/NC_003279.8/chrI/g;s/NC_003280.10/chrII/g;s/NC_003281.10/chrIII/g;s/NC_003282.8/chrIV/g;s/NC_003283.11/chrV/g;s/NC_003284.9/chrX/g'
 

#================================================================#
#                          Report                                #
#================================================================#

NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

#rule report:
#    """
#    Generate a report with the list of datasets + summary of the results.
#    """
# see Scerevisiae report
