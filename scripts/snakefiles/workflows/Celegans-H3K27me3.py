"""This workflow was designed in order to build and test a workflow based on histone marks ChIP-seq analysis.
Data was downloaded fro the GEO platform. It is composed of 2 inputs and 2 ChIP targetting H3K27me3 in C. elegans.
The sequencing was done in single ends.

The workflows should include: 
	- normalization & QC: fastq merging, trimming, quality control...
	- formatting rules: bed, bam, sam, narrowPeak, fasta...
	- alignment: Bowtie, BWA
	- peak-calling: SWEMBL, MACS, SPP, HOMER
	- downstream analyses: IDR, sequence purge, peak length

Usage: snakemake -s scripts/snakefiles/workflows/Celegans-H3K27me3.py -j 12

"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time

#================================================================#
#                    Variables                                   #
#================================================================#

configfile: "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/workflows/Celegans-H3K27me3.json"
workdir: config["dir"]["base"] ## does not work?

# Rules dir
RULES = "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/rules/"

# Raw data
READS = config["dir"]["reads_source"]
RESULTSDIR = config["dir"]["results"]	

CHIP = config["samples"]["chip"].split()
INPUT = config["samples"]["input"].split()
SAMPLES = CHIP + INPUT

GENOME = config["genome"]["genome_version"]

# Blocks to be used /!\ lists vs var
TRIMMING = "sickle_se_q" + config['sickle']['threshold']
ALIGNERS = "bowtie bwa".split()



#================================================================#
#                         Workflow                               #
#================================================================#

# Data conversion & import
for sample in SAMPLES:
	files = os.listdir(READS+sample)
	for f in files:
		raw = f.split(".")[0]
		if not os.path.exists(RESULTSDIR + sample + "/" + sample + ".fastq"):
			os.system("fastq-dump --outdir " + RESULTSDIR + sample + " " + READS + sample + "/" + f)
			os.system("mv " + RESULTSDIR + sample + "/" + raw + ".fastq " + RESULTSDIR + sample + "/" + sample + ".fastq")

# Graphics
GRAPHICS = expand(config["dir"]["results"] + "dag.pdf")

# Data trimming
SICKLE_TRIMMING = expand(config['dir']['results'] + "{samples}/{samples}_" + TRIMMING + ".fastq", samples=SAMPLES)

# Quality control
RAW_QC = expand(config["dir"]["results"] + "{samples}/{samples}_fastqc/", samples=SAMPLES)
TRIMMED_QC = expand(config["dir"]["results"] + "{samples}/{samples}_" + TRIMMING + "_fastqc/", samples=SAMPLES)

# Mapping (can be automatized if several aligners)
#BOWTIE_INDEX = expand()
BOWTIE_MAPPING = expand(config["dir"]["results"] + "{samples}/{samples}_" + TRIMMING + "_bowtie.sam", samples=SAMPLES)

BWA_INDEX = expand(config["dir"]["genome"] + "BWAIndex/{genome}/{genome}.fa.bwt", genome=GENOME)
BWA_MAPPING = expand(config["dir"]["results"] + "{samples}/{samples}_" + TRIMMING + "_bwa.sam", samples=SAMPLES)

# File conversion
SAM_TO_BAM = expand(config["dir"]["results"] + "{sample}/{sample}_" + TRIMMING + "_{aligner}.bam", sample=SAMPLES, aligner=ALIGNERS)
BAM_TO_BED = expand(config["dir"]["results"] + "{sample}/{sample}_" + TRIMMING + "_{aligner}.bed", sample=SAMPLES, aligner=ALIGNERS)

## File conversion
#BED_TO_FASTA = expand(config["dir"]["results"] + '{sample}.fasta', sample=SAMPLES)

## Sequence purge
#PURGED_SEQ = expand(config["dir"]["results"] + "{sample}_purged.fasta", sample=SAMPLES)
##PURGED_SEQ = expand("results/H3K27me3/liver/GSM537698_purged.fasta")

## Oligo analysis
#OLIGO = config['oligo_stats'].split()
#OLIGO_ANALYSIS = expand(config["dir"]["results"] + '{sample}_purged_oligo{oligo}.txt', sample=SAMPLES, oligo=OLIGO)

## Peaks length
#PEAK_LENGTH = expand(config["dir"]["results"] + '{sample}_purged_length.png', sample=SAMPLES)


rule all: 
    """
    Run all the required analyses
    """
    input: GRAPHICS, SICKLE_TRIMMING, RAW_QC, TRIMMED_QC, BOWTIE_MAPPING, BWA_MAPPING, SAM_TO_BAM, BAM_TO_BED
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

#, PEAK_LENGTH, OLIGO_ANALYSIS

# BED_TO_FASTA -> PURGED_SEQUENCE -> (OLIGO_ANALYSIS, PEAK_LENGTH, PEAK_MOTIFS) -> all
# GRAPHICS -> all

#================================================================#
#                         Includes                               #
#================================================================#

#include: RULES + "bed_to_fasta.rules"
include: RULES + "bowtie.rules"
include: RULES + "bwa_se.rules"
include: RULES + "convert_bam_to_bed.rules"
include: RULES + "convert_sam_to_bam.rules"
#include: RULES + "count_oligo.rules"
include: RULES + "fastqc.rules"
include: RULES + "flowcharts.rules"
#include: RULES + "homer.rules"
#include: RULES + "idr.rules"
#include: RULES + "macs14.rules"
#include: RULES + "merge.rules"
#include: RULES + "narrowpeak_to_bed.rules"
#include: RULES + "peak_length.rules"
#include: RULES + "purge_sequence.rules"
include: RULES + "sickle_se.rules"
#include: RULES + "sorted_bam.rules"
#include: RULES + "spp.rules"
#include: RULES + "swembl.rules"
#include: RULES + "trimming.rules"

