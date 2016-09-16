"""Generic work flow for the detection of diferentially expressed
genes from RNA-seq data.

This workflow performs the following treatments: 


Under construction: 
 - convert short read archive files (.sra) into fastq format
 - read quality control with fastQC
 - download the reference genome
 - index the genome for subread-align
 - read mapping with subread-align
 - count the reads per gene with subread featureCounts
 - detection of differentially expressed genes with DESeq2 and/or edgeR

Parameters are specified in a yaml-formatted configuration file.

Usage:
    snakemake -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_se.py \
        --configfile gene-regulation/examples/RNA-seq_GSE71562/RNA-seq_GSE71562.yml -pn

    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/rna-seq_workflow_pe.py \
        --configfile path/to/specific/config_file.yml \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/rna-seq_workflow_pe.py \
        --configfile path/to/specific/config_file.yml \
        --force flowcharts

Sequencing type: 	single end

Author: 		Jeanne Cheneby, Justine Long, Lucie Khamvongsa, Claire Rioualen, Jacques van Helden
Contact: 		Jacques.van-Helden@univ-amu.fr

"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
from datetime import datetime
import pandas as pd

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

if not ("diffexpr" in config["dir"].keys()):
    DEG_DIR = config["dir"]["results"] + "/"
else:
    DEG_DIR = config["dir"]["diffexpr"] + "/"


#================================================================#
#               Snakemake includes
#================================================================#

RULES = os.path.join(GENEREG_LIB, "scripts/snakefiles/rules")
include: os.path.join(RULES, "bam_by_pos.rules")
include: os.path.join(RULES, "bam_stats.rules")
include: os.path.join(RULES, "bam_to_bed.rules")
include: os.path.join(RULES, "bedgraph_to_tdf.rules")
include: os.path.join(RULES, "bowtie.rules")
include: os.path.join(RULES, "bowtie2.rules")
include: os.path.join(RULES, "bowtie2_index.rules")
include: os.path.join(RULES, "bowtie_index.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "cufflinks.rules")
include: os.path.join(RULES, "dot_graph.rules")
include: os.path.join(RULES, "dot_to_image.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "sartools_targetfile.rules")
include: os.path.join(RULES, "genome_coverage_bedgraph.rules")
include: os.path.join(RULES, "get_chrom_sizes.rules")
include: os.path.join(RULES, "gzip.rules")
include: os.path.join(RULES, "index_bam.rules")
include: os.path.join(RULES, "sartools_DESeq2.rules")
include: os.path.join(RULES, "sartools_edgeR.rules")
include: os.path.join(RULES, "sickle.rules")
include: os.path.join(RULES, "sra_to_fastq.rules")
include: os.path.join(RULES, "subread_align.rules")
include: os.path.join(RULES, "subread_featureCounts.rules")
include: os.path.join(RULES, "subread_index.rules")
include: os.path.join(RULES, "tophat.rules")

#================================================================#
#                         Workflow                               #
#================================================================#

## Data import 
IMPORT = expand(FASTQ_DIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS)

## Graphics & reports
GRAPHICS = expand(REPORTS_DIR + "{graph}.png", graph=["dag", "rulegraph"])

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(FASTQ_DIR + "{samples}/{samples}_fastqc/{samples}_fastqc.html", samples=SAMPLE_IDS)
QC = RAW_QC

#----------------------------------------------------------------#
# Trimming
#----------------------------------------------------------------#

if not (("tools" in config.keys()) and ("trimming" in config["tools"].keys())):
    sys.exit("The parameter config['tools']['trimming'] should be specified in the config file. Empty quotes equal to no trimming.")

TRIMMING_TOOLS = config["tools"]["trimming"].split()

TRIMMING = expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}", samples=SAMPLE_IDS, trimmer=TRIMMING_TOOLS)
TRIM = expand("{trimming}.fastq", trimming=TRIMMING)

TRIM_QC = expand(FASTQ_DIR + "{samples}/{samples}_{trimmer}_fastqc/{samples}_{trimmer}_fastqc.html", samples=SAMPLE_IDS, trimmer=TRIMMING_TOOLS)

QC = RAW_QC + TRIM_QC


#----------------------------------------------------------------#
# Alignment
#----------------------------------------------------------------#

if not (("tools" in config.keys()) and ("mapping" in config["tools"].keys())):
    sys.exit("The parameter config['tools']['mapping'] should be specified in the config file.")

MAPPING_TOOLS = config["tools"]["mapping"].split()

INDEX = expand(GENOME_DIR + "{aligner}/" + GENOME_VER + ".fa", aligner=MAPPING_TOOLS)

#MAPPING_TOOLS.append("tophat")

if TRIMMING_TOOLS:
    PREFIX = expand("{trimmer}_{aligner}", aligner=MAPPING_TOOLS, trimmer=TRIMMING_TOOLS)
else:
    PREFIX = expand("{aligner}", aligner=MAPPING_TOOLS)

ALIGNMENT=expand(SAMPLE_DIR + "{samples}/{samples}_{prefix}", samples=SAMPLE_IDS, prefix=PREFIX)

MAPPING = expand("{alignment}.bam", alignment=ALIGNMENT)

SORTED_BAM = expand("{alignment}_sorted_pos.bam", alignment=ALIGNMENT)
SORTED_BAM_BAI = expand("{alignment}_sorted_pos.bam.bai", alignment=ALIGNMENT)
BAM_STATS = expand("{alignment}_bam_stats.txt", alignment=ALIGNMENT)

#----------------------------------------------------------------#
# RNA-seq analysis
#----------------------------------------------------------------#

FEATURE_COUNTS = expand("{alignment}_featureCounts.txt", alignment=ALIGNMENT) ## todo _{feature} , feature=config["subread-featureCounts"]["feature_type"])

INFER_TRANSCRIPTS = expand("{alignment}_cufflinks/transcripts.gtf", alignment=ALIGNMENT)

if not (("tools" in config.keys()) and ("diffexpr" in config["tools"].keys())):
    sys.exit("The parameter config['tools']['diffexpr'] should be specified in the config file.")

DIFFEXPR_TOOLS = config["tools"]["diffexpr"].split()

SARTOOLS_TARGETFILE = expand(DEG_DIR + "{prefix}_SARTools_targetfile.txt", prefix=PREFIX)

DEG = expand(DEG_DIR + "{deg}/{prefix}_{deg}_report.html", prefix=PREFIX, deg=DIFFEXPR_TOOLS)


## ----------------------------------------------------------------
## Visualization & reports
## ----------------------------------------------------------------

## TODO

GENOME_COVERAGE = expand("{alignment}.bedgraph", alignment=ALIGNMENT)

GENOME_COVERAGE_GZ = expand("{alignment}.bedgraph.gz", alignment=ALIGNMENT)
GENOME_COVERAGE_TDF = expand("{alignment}.tdf", alignment=ALIGNMENT)


#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses.
	"""
	input: \
#            IMPORT, \
            QC, \
#            INDEX, \
#            MAPPING, \
#            SORTED_BAM, \
            SORTED_BAM_BAI, \
            BAM_STATS, \
            GENOME_COVERAGE_TDF, \
##            INFER_TRANSCRIPTS, \
#            FEATURE_COUNTS, \
#            SARTOOLS_TARGETFILE, \
            DEG, \
            GRAPHICS
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

##================================================================#
##                         Workflow                               #
##================================================================#

### Data import 
#IMPORT = expand(SAMPLE_DIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS)


### Genome
#GENOME = config["genome"]["version"]
#GENOME_DIR = config["dir"]["genome"] + config["genome"]["version"]

#if not os.path.exists(GENOME_DIR):
#    os.makedirs(GENOME_DIR)

#GENOME_FASTA = expand(GENOME_DIR + "/" + GENOME + ".fa")
#GENOME_ANNOTATIONS = expand(GENOME_DIR + "/" + GENOME + ".{ext}", ext=["gff3", "gtf"])

### Graphics & reports
#GRAPHICS = expand(REPORTS_DIR + "{graph}.png", graph=["dag", "rulegraph"])

##----------------------------------------------------------------#
## Quality control
##----------------------------------------------------------------#

#RAW_QC = expand(SAMPLE_DIR + "{samples}/{samples}_fastqc/{samples}_fastqc.html", samples=SAMPLE_IDS)


##----------------------------------------------------------------#
## Trimming
##----------------------------------------------------------------#

##TRIMMER="sickle-q" + config["sickle"]["threshold"]
##TRIMMING=expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}", samples=SAMPLE_IDS, trimmer=TRIMMER)
##TRIM = expand("{trimming}.fastq", trimming=TRIMMING)

##TRIM_QC = expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}_fastqc/{samples}_{trimmer}_fastqc.html", samples=SAMPLE_IDS, trimmer=TRIMMER)

##QC = RAW_QC + TRIM_QC


##----------------------------------------------------------------#
## Alignment
##----------------------------------------------------------------#


#ALIGNER=["bowtie2", "bowtie", "subread"]
#INDEX = expand(GENOME_DIR + "/{aligner}/" + GENOME + ".fa", aligner=ALIGNER)

#ALIGNER.append("tophat")

#ALIGNMENT=expand(SAMPLE_DIR + "{samples}/{samples}_{aligner}", samples=SAMPLE_IDS, aligner=ALIGNER)
#MAPPING = expand("{alignment}.bam", alignment=ALIGNMENT)

#BAM_STATS = expand("{alignment}_bam_stats.txt", alignment=ALIGNMENT)

##SORTED_BAM = expand("{alignment}_sorted_pos.bam", alignment=ALIGNMENT)

##----------------------------------------------------------------#
## RNA-seq analysis
##----------------------------------------------------------------#


#INFER_TRANSCRIPTS = expand("{alignment}_cufflinks/transcripts.gtf", alignment=ALIGNMENT)


#FEATURE_COUNTS = expand("{alignment}_featureCounts.txt", alignment=ALIGNMENT)

#SARTOOLS_TARGETFILE = expand(DEG_DIR + "{aligner}_SARTools_design.txt", aligner=ALIGNER)

#DEG = expand(DEG_DIR + "{aligner}_{deg}_report.html", aligner=ALIGNER, deg=["DESeq2", "edgeR"])


#GENOME_COVERAGE = expand("{alignment}.bedgraph.gz", alignment=ALIGNMENT)

### ----------------------------------------------------------------
### Visualization
### ----------------------------------------------------------------

##VISU = expand(PEAKS_DIR + "igv_session.xml")

##================================================================#
##                        Rule all                                #
##================================================================#

#rule all: 
#	"""
#	Run all the required analyses.
#	"""
#	input: GRAPHICS, BAM_STATS, FEATURE_COUNTS, INFER_TRANSCRIPTS, GENOME_COVERAGE, SARTOOLS_TARGETFILE, DEG
#	params: qsub=config["qsub"]
#	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

##sed -n FEATURE_COUNTS,  '/digraph/,$p' <rulegraph.dot
