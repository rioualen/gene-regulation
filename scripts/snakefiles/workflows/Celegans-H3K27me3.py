"""This workflow was designed in order to build and test a pipeline based on histone marks ChIP-seq analysis.
Data was downloaded from the GEO platform. It is composed of 2 inputs and 2 ChIP targetting H3K27me3 in C. elegans.

The workflows should include: 
	- normalization & QC: fastq merging, trimming, quality control...
	- formatting rules: bed, bam, sam, narrowPeak, fasta...
	- alignment: Bowtie, BWA
	- peak-calling: SWEMBL, MACS, SPP, HOMER
	- downstream analyses: IDR, sequence purge, peak length

Usage: snakemake -s scripts/snakefiles/workflows/Celegans-H3K27me3.py -j 12 --force flowcharts

Organism: 		Caenorhabditis elegans
Reference genome:	ce10
Sequencing type: 	single end
Data source: 		Gene Expression Omnibus

Author: 		Claire Rioualen
Contact: 		claire.rioualen@inserm.fr
"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time

#================================================================#
#                      Data & directories                        #
#================================================================#

configfile: "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/workflows/Celegans-H3K27me3.json"
workdir: config["dir"]["base"] ## does not work??

# Rules dir
RULES = config["dir"]["rules"]

# Raw data
READS = config["dir"]["reads_source"]
RESULTSDIR = config["dir"]["results"]	

CHIP = config["samples"]["chip"].split()
INPUT = config["samples"]["input"].split()
SAMPLES = CHIP + INPUT

GENOME = config["genome"]["genome_version"]

# /!\ lists vs var
TRIMMING = "sickle_se_q" + config['sickle']['threshold']
ALIGNERS = "bwa".split()


#================================================================#
#                         Workflow                               #
#================================================================#

# Import
IMPORT = expand(RESULTSDIR + "{samples}/{samples}.fastq", samples=SAMPLES)

# Graphics
GRAPHICS = expand(RESULTSDIR + "dag.pdf")

# Data trimming
SICKLE_TRIMMING = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + ".fastq", samples=SAMPLES)

# Quality control
RAW_QC = expand(RESULTSDIR + "{samples}/{samples}_fastqc/", samples=SAMPLES)
TRIMMED_QC = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + "_fastqc/", samples=SAMPLES)

# Mapping (can be automatized if several aligners)
#BOWTIE_INDEX = expand(TODO)
#BOWTIE_MAPPING = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + "_bowtie.sam", samples=SAMPLES)

BWA_INDEX = expand(config["dir"]["genome"] + "BWAIndex/{genome}/{genome}.fa.bwt", genome=GENOME)
BWA_MAPPING = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + "_bwa.sam", samples=SAMPLES)

# File conversion
SAM_TO_BAM = expand(RESULTSDIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.bam", sample=SAMPLES, aligner=ALIGNERS)
BAM_TO_BED = expand(RESULTSDIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.bed", sample=SAMPLES, aligner=ALIGNERS)

# Peak-calling
#HOMER = expand(RESULTSDIR + "{chip}_vs_{inp}/{chip}_vs_{inp}_" + TRIMMING + "_{aligner}_homer.bed", chip=CHIP, inp=INPUT, aligner=ALIGNERS)
#HOMER = expand(RESULTSDIR + "GSM1217459_vs_GSM1217457/GSM1217459_vs_GSM1217457_homer.bed")

MACS2 = expand(RESULTSDIR + "{chip}_vs_{inp}/{chip}_vs_{inp}_" + TRIMMING + "_{aligner}_summits.bed", chip=CHIP, inp=INPUT, aligner=ALIGNERS)

## File conversion
#BED_TO_FASTA = expand(RESULTSDIR + '{sample}.fasta', sample=SAMPLES)

## Sequence purge
#PURGED_SEQ = expand(RESULTSDIR + "{sample}_purged.fasta", sample=SAMPLES)
##PURGED_SEQ = expand("results/H3K27me3/liver/GSM537698_purged.fasta")

## Oligo analysis
#OLIGO = config['oligo_stats'].split()
#OLIGO_ANALYSIS = expand(RESULTSDIR + '{sample}_purged_oligo{oligo}.txt', sample=SAMPLES, oligo=OLIGO)

## Peaks length
#PEAK_LENGTH = expand(RESULTSDIR + '{sample}_purged_length.png', sample=SAMPLES)


rule all: 
    """
    Run all the required analyses
    """
    input: GRAPHICS, RAW_QC, TRIMMED_QC, BWA_MAPPING
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"



# GRAPHICS -> all
# RAW_QC ?
# IMPORT -> SICKLE_TRIMMING -> TRIMMED_QC -> all
# BWA_INDEX -> BWA_MAPPING, BOWTIE_MAPPING -> SAM_TO_BAM -> BAM_TO_BED -> all
# BED_TO_FASTA -> PURGED_SEQUENCE -> (OLIGO_ANALYSIS, PEAK_LENGTH, PEAK_MOTIFS) -> all

rule raw_data:
	"""
	Import & convert sra data to fastq format. 
	/!\ Supposes one file per sample dir /!\
	"""
	input: READS + "{samples}/"
	output: files = RESULTSDIR + "{samples}/{samples}.fastq", \
		dir = RESULTSDIR + "{samples}/"
	params: sample = "{samples}"
	shell: "var=$(ls {input}); fastq-dump --outdir {output.dir} {input}$var; srr=$(basename $var .sra); mv {output.dir}$srr.fastq {output.dir}{params.sample}.fastq"

#================================================================#
#                         Includes                               #
#================================================================#

#include: RULES + "bed_to_fasta.rules"
#include: RULES + "bowtie.rules"
include: RULES + "bwa_index.rules"
include: RULES + "bwa_se.rules"
#include: RULES + "convert_bam_to_bed.rules"
#include: RULES + "convert_sam_to_bam.rules"
#include: RULES + "count_oligo.rules"
include: RULES + "fastqc.rules"
include: RULES + "flowcharts.rules"
#include: RULES + "homer.rules"
#include: RULES + "idr.rules"
#include: RULES + "macs14.rules"
include: RULES + "macs2.rules"
#include: RULES + "merge.rules"
#include: RULES + "narrowpeak_to_bed.rules"
#include: RULES + "peak_length.rules"
#include: RULES + "purge_sequence.rules"
include: RULES + "sickle_se.rules"
#include: RULES + "sorted_bam.rules"
#include: RULES + "spp.rules"
#include: RULES + "swembl.rules"
#include: RULES + "trimming.rules"

