"""Workflow from fastq data to SPP peak caller to IDR calculation...
"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys

#================================================================#
#                     Global variables                           #
#================================================================#

WDIR = "/home/rioualen/Desktop/workspace/fg-chip-seq/"
workdir: WDIR

configfile: "scripts/snakefiles/config_claire.json"

WORKFLOW = config["workflow"]

RESULTSDIR = "results/"

# à automatiser, quand même...
GSM_LIST = "GSM570045 GSM570047".split()
SRR_LIST = "SRR063882 SRR063885".split() 


#CHIP_READ_DIR = config["chip_read_directory"]
CHIP_READ_DIR = "data/chip-seq/reads/"

# Sickle
QUALITY_TYPE = config["sickle_parameters"]["quality_type"]

# Bowtie
BOWTIE_INDEX = config["bowtie_parameters"]["bowtie_index"]
MAX_MISMATCHES = config["bowtie_parameters"]["max_mismatches"]
READS_ISO = config["bowtie_parameters"]["reads_iso"]

ALIGNER = "bowtie".split()
	
#================================================================#
#                         Workflow                               #
#================================================================#



rule all:
	"""Run the workflow on each peak set.
	"""

	input: 
		expand("results/dag.pdf"), \
		expand("results/rulegraph.pdf"), \
		expand("results/{gsm}/{gsm}_merged.fastq", gsm=GSM_LIST), \
		expand("results/{gsm}/{gsm}_merged_fastqc.html", gsm=GSM_LIST), \
		expand("results/{gsm}/{gsm}_trimmed.fastq", gsm=GSM_LIST), \
		expand("results/{gsm}/{gsm}_trimmed_fastqc.html", gsm=GSM_LIST), \
		expand("results/{gsm}/{gsm}_bowtie.sam", gsm=GSM_LIST), \
		expand("results/{gsm}/{gsm}_bwa.sam", gsm=GSM_LIST), \
		expand("results/{gsm}/{gsm}_bowtie.bam", gsm=GSM_LIST)

#		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_purge.fasta", factor=FACTORS, swembl_r=SWEMBL_R), \
#		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_{oligo}.txt", factor=FACTORS, swembl_r=SWEMBL_R, oligo=OLIGO_LENGTH), \
#		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_length.png", factor=FACTORS, swembl_r=SWEMBL_R), \

#rule all:
#    input:expand(CHIP_PEAK_DIR + "/{gsm_treatment}_vs_{gsm_control}_bowtie_swembl_R{swembl_r}/{gsm_treatment}_vs_{gsm_control}_bowtie_swembl_R{swembl_r}.bed", \
#    zip, \
#    gsm_treatment = list_treatment_swembl, \
#    gsm_control =  list_control_swembl, \
#    swembl_r = list_swembl_r)

#================================================================#
#                         Includes                               #
#================================================================#


include: "../rules/flowcharts.rules"
include: "../rules/fastqc.rules"
include: "../rules/fastqc_after.rules"
include: "../rules/merge.rules"
include: "../rules/trimming.rules"
include: "../rules/bowtie.rules"
include: "../rules/bwa.rules"
include: "../rules/convert_sam_to_bam.rules"

#include: "rules/convert_sam_to_bam.rules"
#include: "rules/sorted_bam.rules"
#include: "rules/convert_bam_to_bed.rules"
