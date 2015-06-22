"""ChIP-seq workflow including: 
	- normalization & QC: fastq merging, trimming, quality control...
	- formatting rules: bed, bam, sam, narrowPeak, fasta...
	- alignment: Bowtie, BWA
	- peak-calling: SWEMBL, MACS, SPP, HOMER
	- downstream analyses: IDR, sequence purge, peak length
"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time

#================================================================#
#              Configuration (global viariables)                 #
#================================================================#

#start_time = time.time()

WDIR = "/home/rioualen/Desktop/workspace/fg-chip-seq/"
workdir: WDIR

configfile: "scripts/snakefiles/config_claire.json"

WORKFLOW = config["workflow"]

# Data
CHIP = config["samples"]["chip"].split()
INPUT = config["samples"]["input"].split()
GSM_LIST = CHIP + INPUT

GENOME= config["genome_version"]

CHIP_READ_DIR = config["chip_read_directory"]
RESULTSDIR = "results/"	

# QC
QUALITY_STEP= config["fastQC_steps"].split()

# Trimming
QUALITY_TYPE = config["sickle_parameters"]["quality_type"]

# Alignment
ALIGNER = config["aligners"].split()
BOWTIE_INDEX = config["bowtie_parameters"]["bowtie_index"]
MAX_MISMATCHES = config["bowtie_parameters"]["max_mismatches"]
READS_ISO = config["bowtie_parameters"]["reads_iso"]

# Peak-calling
PEAK_CALLER = config["peak_callers"].split()

# Motif analysis
OLIGO_LENGTH = config["oligo_stats"].split()

#================================================================#
#                         Workflow                               #
#================================================================#

rule all:
	"""Run the workflow on Short Read Run data.
	"""
	input: 
		expand("results/dag.pdf"), \
		expand("results/rulegraph.pdf"), \
		expand("results/{gsm}/{gsm}_{step}_fastqc.html", gsm=GSM_LIST, step=QUALITY_STEP), \
		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_{oligo}.txt", chip=CHIP, inp=INPUT, oligo=OLIGO_LENGTH, aligner=ALIGNER, caller=PEAK_CALLER), \			
		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_length.png", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER), \
#		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}.bed", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER), \
		expand("results/IDR/idr.done")


#================================================================#
#                         Includes                               #
#================================================================#

include: "../rules/bed_to_fasta.rules"
include: "../rules/bowtie.rules"
include: "../rules/bwa.rules"
include: "../rules/convert_bam_to_bed.rules"
include: "../rules/convert_sam_to_bam.rules"
include: "../rules/count_oligo.rules"
include: "../rules/fastqc.rules"
include: "../rules/flowcharts.rules"
include: "../rules/homer.rules"
include: "../rules/idr.rules"
include: "../rules/macs14.rules"
include: "../rules/merge.rules"
include: "../rules/narrowpeak_to_bed.rules"
include: "../rules/peak_length.rules"
include: "../rules/purge_sequence.rules"
include: "../rules/spp.rules"
include: "../rules/swembl.rules"
include: "../rules/trimming.rules"
#include: "rules/sorted_bam.rules"



#print "Program took", time.time() - start_time, "to run"



