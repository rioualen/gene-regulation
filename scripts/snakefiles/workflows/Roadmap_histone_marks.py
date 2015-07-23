"""Workflow designed for the differential analysis of histones marks on sevral human tissues. 
Includes:
	- conversion of bed data
	- purge-sequence (masks redundant sequences, thus saving resources)
	- oligo counting sattistics
	- peak length
	- peak-motifs analysis
	- matrix-quality
	- ...
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

configfile: "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/workflows/histone_marks.json"
workdir: config["dir"]["base"] ## does not work?

# Rules dir
RULES = "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/rules/"

# Raw data, non auto
#MARKS = "H3K27me3 H3K27ac".split()
#TISSUES = "liver heart".split()
#SAMPLES = "GSM537698 GSM910562 GSM1112809 GSM906396".split()

#MARKS = "H3K27me3 H3K27ac".split()
#TISSUES = "liver heart".split()
#SAMPLES = "GSM537698 GSM910562 GSM1112809 GSM906396".split()
#SAMPLES = "H3K27me3/liver/GSM537698 H3K27me3/heart/GSM910562 H3K27ac/liver/GSM1112809 H3K27ac/heart/GSM906396".split()
SAMPLES = "H3K27me3/liver/GSM537698 H3K27me3/heart/GSM910562".split()

# Graph 
GRAPHICS = expand(config["dir"]["results"] + "dag.pdf")

# File conversion
BED_TO_FASTA = expand(config["dir"]["results"] + '{sample}.fasta', sample=SAMPLES)

# Sequence purge
PURGED_SEQ = expand(config["dir"]["results"] + "{sample}_purged.fasta", sample=SAMPLES)
#PURGED_SEQ = expand("results/H3K27me3/liver/GSM537698_purged.fasta")

# Oligo analysis
OLIGO = config['oligo_stats'].split()
OLIGO_ANALYSIS = expand(config["dir"]["results"] + '{sample}_purged_oligo{oligo}.txt', sample=SAMPLES, oligo=OLIGO)

# Peaks length
PEAK_LENGTH = expand(config["dir"]["results"] + '{sample}_purged_length.png', sample=SAMPLES)


rule all: 
    """
    Run all the required analyses
    """
    input: GRAPHICS, BED_TO_FASTA, PURGED_SEQ, PEAK_LENGTH
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

#, PEAK_LENGTH, OLIGO_ANALYSIS

# BED_TO_FASTA -> PURGED_SEQUENCE -> (OLIGO_ANALYSIS, PEAK_LENGTH, PEAK_MOTIFS) -> all
# GRAPHICS -> all

#================================================================#
#                         Includes                               #
#================================================================#

include: RULES + "bed_to_fasta.rules"
#include: RULES + "bowtie.rules"
#include: RULES + "bwa.rules"
#include: RULES + "convert_bam_to_bed.rules"
#include: RULES + "convert_sam_to_bam.rules"
include: RULES + "count_oligo.rules"
#include: RULES + "fastqc.rules"
include: RULES + "flowcharts.rules"
#include: RULES + "homer.rules"
#include: RULES + "idr.rules"
#include: RULES + "macs14.rules"
#include: RULES + "merge.rules"
#include: RULES + "narrowpeak_to_bed.rules"
include: RULES + "peak_length.rules"
include: RULES + "purge_sequence.rules"
#include: RULES + "sickle_se.rules"
#include: RULES + "sorted_bam.rules"
#include: RULES + "spp.rules"
#include: RULES + "swembl.rules"
#include: RULES + "trimming.rules"

