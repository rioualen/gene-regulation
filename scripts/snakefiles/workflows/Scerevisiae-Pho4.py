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

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        --force flowcharts

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

# Beware: verbosity messages are incompatible with the flowcharts
verbosity = int(config["verbosity"])

# Rules dir
RULES = config["dir"]["rules"]

# Raw data
READS = config["dir"]["reads_source"]
RESULTS_DIR = config["dir"]["results"]	
if verbosity >= 1:
    print("\nRESULTS_DIR\t" + RESULTS_DIR)

GENOME = config["genome"]["genome_version"]
if verbosity >= 1:
    print("\nGENOME\t" + GENOME)

# list of suffixes used
#  /!\ can be a *list* of tools?
TRIMMING = "sickle-se-q" + config['sickle']['threshold']
ALIGNER = "bwa"
PEAK_CALLER = "macs2"


#================================================================#
#                         Includes                               #
#================================================================#
if verbosity >= 2:
    print("\nImporting rules")

include: config["dir"]["python_lib"] + "util.py"
include: RULES + "util.rules"
#include: RULES + "assign_samples.rules"
#include: RULES + "bed_to_fasta.rules"
#include: RULES + "bowtie.rules"
include: RULES + "count_reads.rules"
include: RULES + "bwa_index.rules"
include: RULES + "bwa_se.rules"
include: RULES + "clean.rules"
include: RULES + "convert_bam_to_bed.rules"
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


#================================================================#
#                Read sample IDs and design                      #
#================================================================#


## Read sample descriptions
if verbosity >= 1:
    print("\nReading sample descriptions\t" + config["files"]["samples"])
# SAMPLES = read_sample_ids(config["files"]["samples"], verbose=1)
# SAMPLES = TREATMENT + CONTROL ## Read from the sample description file
# Read the sample description file
SAMPLES = read_table(config["files"]["samples"], verbosity=verbosity)
SAMPLE_IDS = SAMPLES.iloc[:,0] ## First column MUST contain the sample ID
SAMPLE_CONDITIONS = SAMPLES.iloc[:,1] ## Second column MUST contain condition for each sample
SAMPLE_EXT_DB = SAMPLES.iloc[:,2] ## External database from which the raw reads will be downloaded
SAMPLE_EXT_ID = SAMPLES.iloc[:,3] ## Sample ID in the external database
SAMPLE_DESCR = SAMPLES.iloc[:,7] ## First column MUST contain the sample ID
if verbosity >= 2:
    print("\nSAMPLE_IDS\t" + ";".join(SAMPLE_IDS))

## Read the list of sample IDs from the sample description file
if verbosity >= 1:
    print("\nReading experimental design\t" + config["files"]["design"])
# TREATMENT = config["samples"]["chip"].split()
# CONTROL = config["samples"]["input"].split()
DESIGN = read_table(config["files"]["design"], verbosity=verbosity)
TREATMENT = DESIGN.iloc[:,0]
CONTROL = DESIGN.iloc[:,1]
#TREATMENT, CONTROL = read_chipseq_design(config["files"]["design"], test_column=1, input_column=2, verbose=2)
if verbosity >= 2:
    print("\nTREATMENT\t" + ";".join(TREATMENT))
    print("\nCONTROL\t" + ";".join(CONTROL))

#================================================================#
#                         Workflow                               #
#================================================================#

# ## Data import & merging.
# ##
# ## Tricky python code to prepare the data, but this should be
# ## commented to avoid executing it at each run of the workflow. Should
# ## be converted to a rule.
# if not os.path.exists(RESULTS_DIR):
# 	os.makedirs(RESULTS_DIR)
# for sample in SAMPLES:
# 	indir = READS + sample + "/"
# 	outdir = RESULTS_DIR + sample + "/"
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
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")

# Data trimming
TRIMMED_READS_SICKLE = expand(RESULTS_DIR + "{samples}/{samples}_" + TRIMMING + ".fastq", samples=SAMPLE_IDS)
if verbosity >= 3:
    print("\nTRIMMED_READS_SICKLE\n\t" + "\n\t".join(TRIMMED_READS_SICKLE))

# Quality control
RAW_QC = expand(RESULTS_DIR + "{samples}/{samples}_fastqc/", samples=SAMPLE_IDS)
if verbosity >= 3:
    print("\nRAW_QC\n\t" + "\n\t".join(RAW_QC))
RAW_READNB = expand(RESULTS_DIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)
TRIMMED_QC = expand(RESULTS_DIR + "{samples}/{samples}_" + TRIMMING + "_fastqc/", samples=SAMPLE_IDS)
if verbosity >= 3:
    print("\nTRIMMED_QC\n\t" + "\n\t".join(TRIMMED_QC))

# Mapping with BWA
BWA_INDEX = expand(config["dir"]["genome"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
MAPPED_READS_BWA = expand(RESULTS_DIR + "{samples}/{samples}_" + TRIMMING + "_{aligner}.sam", samples=SAMPLE_IDS, aligner=ALIGNER)
if verbosity >= 3:
    print("\nMAPPED_READS_BWA\n\t" + "\n\t".join(MAPPED_READS_BWA))

# Sorted and converted reads (bam, bed)
SORTED_MAPPED_READS_BWA = expand(RESULTS_DIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}_sorted_pos.bam", sample=SAMPLE_IDS, aligner=ALIGNER)
BAM_READNB = expand(RESULTS_DIR + "{samples}/{samples}_" + TRIMMING + "_{aligner}_sorted_pos_bam_readnb.txt", samples=SAMPLE_IDS, aligner=ALIGNER)
SORTED_READS_BED = expand(RESULTS_DIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}_sorted_pos.bed", sample=SAMPLE_IDS, aligner=ALIGNER)
BED_READNB = expand(RESULTS_DIR + "{samples}/{samples}_" + TRIMMING + "_{aligner}_sorted_pos_bed_nb.txt", samples=SAMPLE_IDS, aligner=ALIGNER)
if verbosity >= 3: 
    print("\nSORTED_READS_BED\n\t" + "\n\t".join(SORTED_READS_BED))
#CONVERTED_BED = expand(RESULTS_DIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.converted.bed", sample=SAMPLE_IDS, aligner=ALIGNER)

# Peak-calling
# ! In case of use of several peak-callers, beware of specific prefixes or such...
#MACS2 = expand(RESULTS_DIR + "{treat}_vs_{ctrl}/{treat}_vs_{ctrl}_{trimming}_{aligner}_{caller}_peaks.narrowPeak", treat="GSM121459", ctrl="GSM1217457", trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
MACS2 = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/{treat}_vs_{ctrl}_{{trimming}}_{{aligner}}_{{caller}}_summits.bed", 
               zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
if verbosity >= 3: 
    print("\nMACS2\n\t" + "\n\t".join(MACS2))

# File conversion
#NARROWPEAK_TO_BED = expand(RESULTS_DIR + "{treat}_vs_{ctrl}/{treat}_vs_{ctrl}_{trimming}_{aligner}_{caller}_peaks.bed", treat=TREATMENT, ctrl=CONTROL, trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
#BED_TO_FASTA = expand(RESULTS_DIR + "{treat}_vs_{ctrl}/{treat}_vs_{ctrl}_{trimming}_{aligner}_{caller}_peaks.fasta", treat=TREATMENT, ctrl=CONTROL, trimming=TRIMMING, aligner=ALIGNER, caller=PEAK_CALLER)
BED_TO_FASTA = expand(RESULTS_DIR + "{sample}/{sample}_" + TRIMMING + "_{aligner}.calling.fasta", sample=SAMPLE_IDS, trimming=TRIMMING, aligner=ALIGNER)

## Sequence purge
#PURGED_SEQ = expand(RESULTS_DIR + "{sample}_purged.fasta", sample=SAMPLE_IDS)
##PURGED_SEQ = expand("results/H3K27me3/liver/GSM537698_purged.fasta")

## Oligo analysis
#OLIGO = config['oligo_stats'].split()
#OLIGO_ANALYSIS = expand(RESULTS_DIR + '{sample}_purged_oligo{oligo}.txt', sample=SAMPLE_IDS, oligo=OLIGO)

## Peaks length
#PEAK_LENGTH = expand(RESULTS_DIR + '{sample}_purged_length.png', sample=SAMPLE_IDS)

ruleorder: macs2 > bam_to_bed > sam2bam
ruleorder: count_reads_bam > sam2bam

CLEANING = expand(RESULTS_DIR + "cleaning.done")

rule all: 
    """
    Run all the required analyses
    """
#    input: GRAPHICS, RAW_READNB, RAW_QC, TRIMMED_QC, SORTED_MAPPED_READS_BWA, SORTED_READS_BED
#    input: RAW_QC, TRIMMED_QC, MAPPED_READS_BWA, BAM_READNB, BED_READNB
    input: MACS2 ### NOT WORKING YET
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"


#================================================================#
#                       Local rules                              #
#================================================================#

