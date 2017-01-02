
"""Generic work< flow for the analysis of ChIP-seq data for the binding
of transcription factors.


This workflow performs the following treatments: 

 - read mapping
 - peak-calling with alternate peak-calling programs
( - motif discovery)

The details are specified in a yaml-formatted configuration file.

Usage:

1. Run in command line mode

    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/factor_workflow.py 
        --configfile path/to/specific/config_file.yml \
        [targets]

2. Send tasks to qsub job scheduler

    snakemake -p -c "qsub {params.qsub}" -j 12 \
        -s gene-regulation/scripts/snakefiles/workflows/factor_workflow.py \
        --configfile path/to/specific/config_file.yml \
        [targets]

3. Generate flowcharts:

    snakemake -p -s scripts/snakefiles/workflows/factor_workflow.py \
        --configfile path/to/specific/config_file.yml \
        --force flowcharts

Sequencing type: 	single end

Author: 		Claire Rioualen, Jacques van Helden
Contact: 		claire.rioualen@inserm.fr

"""



#================================================================#
#                       Python Imports 
#================================================================#

from snakemake.utils import R
import os
import sys
import datetime
import re
import pandas as pd

# Define verbosity
if not ("verbosity" in config.keys()):
    config["verbosity"] = 0
verbosity = int(config["verbosity"])


GENEREG_LIB = os.path.abspath(config["dir"]["gene_regulation"])

# Python includes
PYTHON = os.path.join(GENEREG_LIB, "scripts/python_lib")
include: os.path.join(PYTHON, "util.py")

#================================================================#
#                      Global variables
#================================================================#

# Raw data
if not (("dir" in config.keys()) and ("reads_source" in config["dir"].keys())):
    sys.exit("The parameter config['dir']['reads_source'] should be specified in the config file.")
else:
    READS = config["dir"]["reads_source"]

# Samples
SAMPLES = read_table(config["metadata"]["samples"])
SAMPLE_IDS = SAMPLES.iloc[:,0]
SAMPLE_CONDITIONS = read_table(config["metadata"]["samples"])['condition']

# Genome & annotations
GENOME_VER = config["genome"]["version"]
GENOME_DIR = config["dir"]["genome"] + "/"
GENOME_FASTA = GENOME_DIR + config["genome"]["fasta_file"]
GENOME_GFF3 = GENOME_DIR + config["genome"]["gff3_file"]
GENOME_GTF = GENOME_DIR + config["genome"]["gtf_file"]

# Design
DESIGN = read_table(config["metadata"]["design"])
TREATMENT = DESIGN['treatment']
CONTROL = DESIGN['control']

## Data & results dir
if not (("dir" in config.keys()) and ("reads_source" in config["dir"].keys())):
    sys.exit("The parameter config['dir']['reads_source'] should be specified in the config file.")

READS = config["dir"]["reads_source"]
if not os.path.exists(READS):
    os.makedirs(READS)

if not ("results" in config["dir"].keys()):
    sys.exit("The parameter config['dir']['results'] should be specified in the config file.")
else:
    RESULTS_DIR = config["dir"]["results"] + "/"

if not ("fastq" in config["dir"].keys()):
    sys.exit("The parameter config['dir']['fastq'] should be specified in the config file.")
else:
    FASTQ_DIR = config["dir"]["fastq"] + "/"

if not ("samples" in config["dir"].keys()):
    SAMPLE_DIR = config["dir"]["results"] + "/"
else:
    SAMPLE_DIR = config["dir"]["samples"] + "/"

if not ("reports" in config["dir"].keys()):
    REPORTS_DIR = config["dir"]["results"] + "/"
else:
    REPORTS_DIR = config["dir"]["reports"] + "/"

if not ("peaks" in config["dir"].keys()):
    PEAKS_DIR = config["dir"]["results"] + "/"
else:
    PEAKS_DIR = config["dir"]["peaks"] + "/"


#================================================================#
#               Snakemake includes
#================================================================#

RULES = os.path.join(GENEREG_LIB, "scripts/snakefiles/rules")

include: os.path.join(RULES, "bam_by_pos.rules")
include: os.path.join(RULES, "bam_to_bed.rules")
include: os.path.join(RULES, "bam_stats.rules")
include: os.path.join(RULES, "bowtie_index.rules")
include: os.path.join(RULES, "bowtie.rules")
include: os.path.join(RULES, "bowtie2_index.rules")
include: os.path.join(RULES, "bowtie2.rules")
include: os.path.join(RULES, "bPeaks.rules")
include: os.path.join(RULES, "bwa_index.rules")
include: os.path.join(RULES, "bwa.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "dot_graph.rules")
include: os.path.join(RULES, "dot_to_image.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "genome_coverage_bedgraph.rules")
include: os.path.join(RULES, "getfasta.rules")
include: os.path.join(RULES, "get_chrom_sizes.rules")
include: os.path.join(RULES, "gzip.rules")
include: os.path.join(RULES, "homer.rules")
include: os.path.join(RULES, "igv_session.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "macs14.rules")
include: os.path.join(RULES, "peak_motifs.rules")
include: os.path.join(RULES, "sickle.rules")
include: os.path.join(RULES, "spp.rules")
include: os.path.join(RULES, "subread_index.rules")
include: os.path.join(RULES, "subread.rules")
include: os.path.join(RULES, "swembl.rules")
include: os.path.join(RULES, "sra_to_fastq.rules")


#================================================================#
#                      Data & wildcards                          #
#================================================================#

# Raw data
READS = config["dir"]["reads_source"]

# Samples
SAMPLES = read_table(config["metadata"]["samples"])
SAMPLE_IDS = SAMPLES.iloc[:,0]


## Design
DESIGN = read_table(config["metadata"]["design"])
TREATMENT = DESIGN['treatment']
CONTROL = DESIGN['control']

## Data & results dir

if not (("dir" in config.keys()) and ("reads_source" in config["dir"].keys())):
    sys.exit("The parameter config['dir']['reads_source'] should be specified in the config file.")

READS = config["dir"]["reads_source"]
if not os.path.exists(READS):
    os.makedirs(READS)

if not ("results" in config["dir"].keys()):
    sys.exit("The parameter config['dir']['results'] should be specified in the config file.")

RESULTS_DIR = config["dir"]["results"]
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

if not ("samples" in config["dir"].keys()):
    SAMPLE_DIR = config["dir"]["results"]
else:
    SAMPLE_DIR = config["dir"]["samples"]
if not os.path.exists(SAMPLE_DIR):
    os.makedirs(SAMPLE_DIR)

if not ("peaks" in config["dir"].keys()):
    PEAKS_DIR = config["dir"]["results"]
else:
    PEAKS_DIR = config["dir"]["peaks"]
if not os.path.exists(PEAKS_DIR):
    os.makedirs(PEAKS_DIR)

if not ("reports" in config["dir"].keys()):
    REPORTS_DIR = config["dir"]["results"]
else:
    REPORTS_DIR = config["dir"]["reports"]
if not os.path.exists(REPORTS_DIR):
    os.makedirs(REPORTS_DIR)



#================================================================#
#                         Workflow                               #
#================================================================#

## Data import 

IMPORT = expand(SAMPLE_DIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS)

# Genome
GENOME = config["genome"]["version"]
GENOME_DIR = config["dir"]["genome"] + config["genome"]["version"]

if not os.path.exists(GENOME_DIR):
    os.makedirs(GENOME_DIR)

GENOME_FASTA = expand(GENOME_DIR + "/" + GENOME + ".fa")
GENOME_ANNOTATIONS = expand(GENOME_DIR + "/" + GENOME + ".gff3")

#----------------------------------------------------------------
# Quality control
#----------------------------------------------------------------

RAW_QC = expand(SAMPLE_DIR + "{samples}/{samples}_fastqc/{samples}_fastqc.html", samples=SAMPLE_IDS)


#----------------------------------------------------------------
# Trimming
#----------------------------------------------------------------

TRIMMER="sickle-se-q" + config["sickle"]["threshold"]#trimmomatic
TRIMMING=expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}", samples=SAMPLE_IDS, trimmer=TRIMMER)
TRIM = expand(SAMPLE_DIR + "{trimming}.fastq", trimming=TRIMMING)

TRIM_QC = expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}_fastqc/{samples}_{trimmer}_fastqc.html", samples=SAMPLE_IDS, trimmer=TRIMMER)

QC = RAW_QC + TRIM_QC


#----------------------------------------------------------------
# Alignment
#----------------------------------------------------------------


ALIGNER=["bowtie"]# bwa, bowtie2, subread
ALIGNMENT=expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}_{aligner}", samples=SAMPLE_IDS, aligner=ALIGNER, trimmer=TRIMMER)

INDEX = expand(GENOME_DIR + "/{aligner}/" + GENOME + ".fa", aligner=ALIGNER)

MAPPING = expand("{alignment}.sam", alignment=ALIGNMENT)

BAM_STATS = expand("{alignment}_bam_stats.txt", alignment=ALIGNMENT)

GENOME_COVERAGE = expand("{alignment}.bedgraph", alignment=ALIGNMENT)
GENOME_COVERAGE_GZ = expand("{alignment}.bedgraph.gz", alignment=ALIGNMENT)


# Sort mapped reads

#SORTED_READS_BED = expand("{alignment}_sorted_pos.bed", alignment=ALIGNMENT)


# ----------------------------------------------------------------
# Peak-calling
# ----------------------------------------------------------------

PEAKCALLER=[
    "homer-fdr" + config["homer"]["fdr"],
    "macs2-qval" + config["macs2"]["qval"], 
    "swembl-R" + config["swembl"]["R"],
    "macs14-pval" + config["macs14"]["pval"],
    "spp-fdr" + config["spp"]["fdr"],
    "bPeaks-log" + config["bPeaks"]["log2FC"],
]

PEAKCALLING=expand(expand(PEAKS_DIR + "{treat}_vs_{control}/{{peakcaller}}/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER, trimmer=TRIMMER)

PEAKS = expand("{peakcalling}.bed", peakcalling=PEAKCALLING)

# ----------------------------------------------------------------
# Peak annotation
# ----------------------------------------------------------------

GENE_ANNOT = ["intersect", "window"]#"closest"
PEAKS_TO_GENES = expand("{peakcalling}_{gene_annotation}_annot.bed", peakcalling=PEAKCALLING, gene_annotation=GENE_ANNOT)

MOTIFS=expand(expand("{treat}_vs_{control}/{{peakcaller}}/peak-motifs/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}_peak-motifs_synthesis", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER, trimmer=TRIMMER)


GET_FASTA = expand(PEAKS_DIR + "{peakcalling}.fasta", peakcalling=PEAKCALLING)
PEAK_MOTIFS = expand(PEAKS_DIR + "{motifs}.html", motifs=MOTIFS)

# ----------------------------------------------------------------
# Reports
# ----------------------------------------------------------------


GRAPHICS = expand(REPORTS_DIR + "{graph}.{ext}", graph=["dag", "rulegraph"], ext=["png", "pdf"])


## Following not yet properly implemented
#IGV = expand(REPORTS_DIR + "igv_session.xml")
#BED_INTER = expand(REPORTS_DIR + "multiinter.tab")

#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses.
	"""
	input: BAM_STATS, GENOME_COVERAGE_GZ, GRAPHICS, QC, PEAKS_TO_GENES, PEAKS#, PEAK_MOTIFS
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"
