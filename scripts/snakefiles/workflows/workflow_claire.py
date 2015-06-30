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
#              Configuration (global variables)                  #
#================================================================#

#start_time = time.time()

WDIR = "/home/rioualen/Desktop/workspace/fg-chip-seq/"
workdir: WDIR

configfile: "scripts/snakefiles/config_claire.json"

# Data
CHIP = config["samples"]["chip"].split()
INPUT = config["samples"]["input"].split()
GSM_LIST = CHIP + INPUT

READS = config["chip_read_directory"]
RESULTSDIR = config["results_directory"]	
LOGS = config["log_directory"]

# QC
QUALITY_STEP= config["fastQC_steps"].split()

# Alignment
#ALIGNER = config["aligners"].split()
#BOWTIE_INDEX = config["bowtie_parameters"]["bowtie_index"]
#MAX_MISMATCHES = config["bowtie_parameters"]["max_mismatches"]
#READS_ISO = config["bowtie_parameters"]["reads_iso"]

## Peak-calling
#PEAK_CALLER = config["peak_callers"].split()

## Motif analysis
#OLIGO_LENGTH = config["oligo_stats"].split()
## Peak-motifs params TODO
#PEAK_MOTIFS_TASKS = config["peak-motifs_parameters"]["tasks"]
#MAX_SEQ_LENGTH = config["peak-motifs_parameters"]["max_seq_length"]
#NMOTIFS = config["peak-motifs_parameters"]["nmotifs"]
#MINOL = config["peak-motifs_parameters"]["minol"]
#MAXOL = config["peak-motifs_parameters"]["maxol"]
#MOTIF_DB = config["peak-motifs_parameters"]["motif_db"]

#================================================================#
#                         Workflow                               #
#================================================================#

rule all:
	"""Run the workflow on Short Read Run data.
	"""
	input: 
		expand(RESULTSDIR + "dag.pdf"), \
		expand(RESULTSDIR + "rule.pdf"), \
		expand(RESULTSDIR + "{gsm}/{gsm}_trimmed.fastq", gsm=GSM_LIST), \
		expand(RESULTSDIR + "{gsm}/{gsm}_{step}_fastqc.html", gsm=GSM_LIST, step=QUALITY_STEP)
#		expand(RESULTSDIR + "{gsm}/{gsm}_{trimmer}.fastq", gsm=GSM_LIST, trimmer="sickle"), \
#		expand(RESULTSDIR + "{gsm}/{gsm}_{trimmer}_fastqc.html", gsm=GSM_LIST, trimmer="sickle")
##		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_{oligo}.txt", chip=CHIP, inp=INPUT, oligo=OLIGO_LENGTH, aligner=ALIGNER, caller=PEAK_CALLER), \			
#		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_length.png", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER)
#		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}.bed", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER)
#		expand("results/IDR/idr.done")

#================================================================#
#                        Data import                             #
#================================================================#

rule merge_reads:
    """Merge the raw short read files (SRR) belonging to the same sample (GSM)"""
    output: RESULTSDIR + "{gsm}/{gsm}_import.fastq"
    params: qsub = config["qsub"] + " -e " + LOGS + "/{gsm}/{gsm}_merged_qsub.err -o " + LOGS + "/{gsm}/{gsm}_merged_qsub.out"
    log: LOGS + "{gsm}/{gsm}_merged.log"
    benchmark: LOGS + "{gsm}/{gsm}_merged_benchmark.json"
    message: "Merging short reads for sample {input}"
    shell: "ls -1 {READS}/{wildcards.gsm}/SRR*.fastq | xargs cat > {output};"

#================================================================#
#                         Includes                               #
#================================================================#

#include: "../rules/bed_to_fasta.rules"
#include: "../rules/bowtie.rules"
#include: "../rules/bwa.rules"
#include: "../rules/convert_bam_to_bed.rules"
#include: "../rules/convert_sam_to_bam.rules"
#include: "../rules/count_oligo.rules"
include: "../rules/fastqc.rules"
include: "../rules/flowcharts.rules"
#include: "../rules/homer.rules"
#include: "../rules/idr.rules"
#include: "../rules/macs14.rules"
#include: "../rules/merge.rules"
#include: "../rules/narrowpeak_to_bed.rules"
#include: "../rules/peak_length.rules"
#include: "../rules/purge_sequence.rules"
include: "../rules/sickle_se.rules"
#include: "../rules/sorted_bam.rules"
#include: "../rules/spp.rules"
#include: "../rules/swembl.rules"
#include: "../rules/trimming.rules"




#print "Program took", time.time() - start_time, "to run"



