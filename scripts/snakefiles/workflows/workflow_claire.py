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
#                    Variables                                   #
#================================================================#

configfile: "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/config_claire.json"
workdir: config["dir"]["base"] ## does not work?

# Data
CHIP = config["samples"]["chip"].split()
INPUT = config["samples"]["input"].split()
SAMPLES = CHIP + INPUT

READS = config["dir"]["reads"]
RESULTSDIR = config["dir"]["results"]	
LOGS = config["dir"]["logs"]

# Rules dir
RULES = "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/rules/"

# QC
QUALITY_STEP= config["fastqc"]["steps"].split()

#ALN
ALIGNER = config['aligners'].split()


#================================================================#
#                         Workflow                               #
#================================================================#

rule all:
	"""Run the workflow on Short Read Run data.
	"""
	input: 
		expand(RESULTSDIR + "dag.pdf"), \
		expand(RESULTSDIR + "rule.pdf"), \
		expand(RESULTSDIR + "{samples}/{samples}_sickle.fastq", samples=SAMPLES), \
		expand(RESULTSDIR + "{samples}/{samples}_{step}_fastqc.html", samples=SAMPLES, step=QUALITY_STEP), \
		expand(RESULTSDIR + "{samples}/{samples}_{aligner}.sam", samples=SAMPLES, aligner=ALIGNER)
#		expand(RESULTSDIR + "{samples}/{samples}_{trimmer}.fastq", gsm=GSM_LIST, trimmer="sickle"), \
#		expand(RESULTSDIR + "{samples}/{samples}_{trimmer}_fastqc.html", gsm=GSM_LIST, trimmer="sickle")
##		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_{oligo}.txt", chip=CHIP, inp=INPUT, oligo=OLIGO_LENGTH, aligner=ALIGNER, caller=PEAK_CALLER), \			
#		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_length.png", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER)
#		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}.bed", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER)
#		expand("results/IDR/idr.done")

#================================================================#
#                        Data import                             #
#================================================================#

rule merge_reads:
    """Merge the raw short read files (SRR) belonging to the same sample (GSM) before allowing the trimming."""
    output: RESULTSDIR + "{samples}/{samples}_import.fastq"
    params: qsub = config["qsub"] + " -e " + LOGS + "/{samples}/{samples}_merged_qsub.err -o " + LOGS + "/{samples}/{samples}_merged_qsub.out"
    log: LOGS + "{samples}/{samples}_merged.log"
    benchmark: LOGS + "{samples}/{samples}_merged_benchmark.json"
    message: "Merging short reads for sample {input}"
    shell: "ls -1 {READS}/{wildcards.gsm}/SRR*.fastq | xargs cat > {output};"

#================================================================#
#                         Includes                               #
#================================================================#

#include: RULES + "bed_to_fasta.rules"
include: RULES + "bowtie.rules"
include: RULES + "bwa.rules"
#include: RULES + "convert_bam_to_bed.rules"
#include: RULES + "convert_sam_to_bam.rules"
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

