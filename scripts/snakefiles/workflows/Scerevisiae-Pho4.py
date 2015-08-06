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

Usage: snakemake -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py -p  -c "qsub {params.qsub}" -j 12 --force flowcharts

Organism: 		Saccharomyces cerevisiae
Reference genome:	
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

configfile: "scripts/snakefiles/workflows/Scerevisiae-Pho4.json"
#workdir: config["dir"]["base"] ## does not work??


# Rules dir
RULES = config["dir"]["rules"]

# Raw data
READS = config["dir"]["reads_source"]
RESULTSDIR = config["dir"]["results"]	

# CHIP = config["samples"]["chip"].split()
# INPUT = config["samples"]["input"].split()
# SAMPLES = CHIP + INPUT ## Read from the sample description file

GENOME = config["genome"]["genome_version"]

# list of suffixes used
#  /!\ can be a *list* of tools?
TRIMMING = "sickle-se-q" + config['sickle']['threshold']
ALIGNER = "bwa"
PEAK_CALLER = "macs2"


#================================================================#
#                         Includes                               #
#================================================================#

include: config["dir"]["python_lib"] + "util.py"
#include: RULES + "assign_samples.rules"
#include: RULES + "bed_to_fasta.rules"
#include: RULES + "bowtie.rules"
include: RULES + "bwa_index.rules"
include: RULES + "bwa_se.rules"
include: RULES + "clean.rules"
include: RULES + "convert_bam_to_bed.rules"
include: RULES + "convert_sam_to_bam.rules"
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


#================================================================#
#                Read sample IDs and design                      #
#================================================================#
## Read the list of sample IDs from the sample description file
SAMPLES = read_sample_ids(config["files"]["samples"], verbose=1)
CHIP, INPUT = read_chipseq_design(config["files"]["design"], test_column=1, input_column=2, verbose=2)

#================================================================#
#                         Workflow                               #
#================================================================#

# ## Data import & merging.
# ##
# ## Tricky python code to prepare the data, but this should be
# ## commented to avoid executing it at each run of the workflow. Should
# ## be converted to a rule.
# if not os.path.exists(RESULTSDIR):
# 	os.makedirs(RESULTSDIR)
# for sample in SAMPLES:
# 	indir = READS + sample + "/"
# 	outdir = RESULTSDIR + sample + "/"
# 	if not os.path.exists(outdir):
# 		os.makedirs(outdir)
# 	sra_files = os.listdir(indir)
# 	for sra in sra_files:
# 		os.system("fastq-dump --outdir " + outdir + " " + indir + sra)
# 	if len(sra_files) > 1:
# 		os.system("ls -1 " + outdir + "*.fastq | xargs cat > " + outdir + sample + ".fastq")
# 		os.system("find " + outdir + " -type f -name SRR* | xargs rm")
# 	else:
# 		fastq_files = os.listdir(outdir)
# 		cmd = "ls -1 " + outdir + fastq_files[0] + " | xargs mv " + outdir + sample + ".fastq"
# 		print("Running command: "+ cmd)
# 		os.system(cmd)
	

# Graphics
GRAPHICS = expand(RESULTSDIR + "dag.pdf")

# Data trimming
SICKLE_TRIMMING = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + ".fastq", samples=SAMPLES)

# Quality control
RAW_QC = expand(RESULTSDIR + "{samples}/{samples}_fastqc/", samples=SAMPLES)
TRIMMED_QC = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + "_fastqc/", samples=SAMPLES)

# Mapping
BWA_INDEX = expand(config["dir"]["genome"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
BWA_MAPPING = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + "_{aligner}.sam", samples=SAMPLES, aligner=ALIGNER)

# File conversion
SAM_TO_BAM = expand(RESULTSDIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.bam", sample=SAMPLES, aligner=ALIGNER)
#BAM_COUNTS = expand(RESULTSDIR + "{samples}/{samples}_" + TRIMMING + "_{aligner}_bam_reandb.txt", samples=SAMPLES, aligner=ALIGNER)
BAM_TO_BED = expand(RESULTSDIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.bed", sample=SAMPLES, aligner=ALIGNER)
print("BAM_TO_BED\n\t" + "\n\t".join(BAM_TO_BED))
#CONVERTED_BED = expand(RESULTSDIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.converted.bed", sample=SAMPLES, aligner=ALIGNER)

# Peak-calling
# ! In case of use of several peak-callers, beware of specific prefixes or such...
#MACS2 = expand(RESULTSDIR + "{chip}_vs_{inp}/{chip}_vs_{inp}_{trimming}_{aligner}_{caller}_peaks.narrowPeak", chip="GSM121459", inp="GSM1217457", trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
MACS2 = expand(expand(RESULTSDIR + "{chip}_vs_{inp}/{chip}_vs_{inp}_{{trimming}}_{{aligner}}_{{caller}}._summits.bed", 
               zip, chip=CHIP, inp=INPUT), trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
print("MACS2\n\t" + "\n\t".join(MACS2))

# File conversion
#NARROWPEAK_TO_BED = expand(RESULTSDIR + "{chip}_vs_{inp}/{chip}_vs_{inp}_{trimming}_{aligner}_{caller}_peaks.bed", chip=CHIP, inp=INPUT, trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
#BED_TO_FASTA = expand(RESULTSDIR + "{chip}_vs_{inp}/{chip}_vs_{inp}_{trimming}_{aligner}_{caller}_peaks.fasta", chip=CHIP, inp=INPUT, trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
BED_TO_FASTA = expand(RESULTSDIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.calling.fasta", sample=SAMPLES, trimming=TRIMMING, aligner=ALIGNER)

## Sequence purge
#PURGED_SEQ = expand(RESULTSDIR + "{sample}_purged.fasta", sample=SAMPLES)
##PURGED_SEQ = expand("results/H3K27me3/liver/GSM537698_purged.fasta")

## Oligo analysis
#OLIGO = config['oligo_stats'].split()
#OLIGO_ANALYSIS = expand(RESULTSDIR + '{sample}_purged_oligo{oligo}.txt', sample=SAMPLES, oligo=OLIGO)

## Peaks length
#PEAK_LENGTH = expand(RESULTSDIR + '{sample}_purged_length.png', sample=SAMPLES)

ruleorder: macs2 > bam_to_bed

CLEANING = expand(RESULTSDIR + "cleaning.done")

rule all: 
    """
    Run all the required analyses
    """
#    input: GRAPHICS, RAW_QC, TRIMMED_QC, BAM_TO_BED
    input: MACS2
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"


#================================================================#
#                       Local rules                              #
#================================================================#

