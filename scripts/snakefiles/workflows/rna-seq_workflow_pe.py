"""Generic work flow for the detection of diferentially expressed
genes from RNA-seq data using paired-ends sequencing type.

This workflow performs the following treatments: 

Under construction: 
 - convert short read archive files (.sra) into fastq format
 - read quality control with fastQC
 - download the reference genome
 - index the genome for subread-align
 - read mapping with subread-align (other possible aligners include: bowtie, bowtie2, tophat)
 - count the reads per gene with subread featureCounts
 - detection of differentially expressed genes with DESeq2 and/or edgeR using the SARTools package

Parameters are specified in a yaml-formatted configuration file.

Usage:
    snakemake -p -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_pe.py --configfile gene-regulation/examples/RNA-seq_GSE41190/RNA-seq_GSE41190.yml

    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_pe.py \
        --configfile gene-regulation/examples/RNA-seq_GSE41190/RNA-seq_GSE41190.yml \
        [targets]

Sequencing type: 	paired end

Author: 		Jeanne Cheneby, Justine Long, Lucie Khamvongsa, Claire Rioualen, Jacques van Helden
Contact: 		Jacques.van-Helden@univ-amu.fr

"""

#================================================================#
#                        Imports & includes
#================================================================#

from snakemake.utils import R
import os
import sys
import datetime
import pandas as pd


#if not ("base" in config["dir"].keys()):
#    sys.exit("The parameter 'dir base' should be specified in the config file.")

GENEREG_LIB = os.path.abspath(config["dir"]["gene_regulation"])

# Snakemake includes
RULES = os.path.join(GENEREG_LIB, "scripts/snakefiles/rules")
include: os.path.join(RULES, "annotation_download.rules")
include: os.path.join(RULES, "bam_by_pos.rules")
include: os.path.join(RULES, "bam_stats.rules")
include: os.path.join(RULES, "bam_to_bed.rules")
include: os.path.join(RULES, "bowtie.rules")
include: os.path.join(RULES, "bowtie_index.rules")
include: os.path.join(RULES, "bowtie2.rules")
include: os.path.join(RULES, "bowtie2_index.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "cufflinks.rules")
include: os.path.join(RULES, "dot_graph.rules")
include: os.path.join(RULES, "dot_to_image.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "genome_download.rules")
include: os.path.join(RULES, "get_chrom_sizes.rules")
include: os.path.join(RULES, "sartools_DESeq2.rules")
include: os.path.join(RULES, "sartools_edgeR.rules")
include: os.path.join(RULES, "sartools_targetfile.rules")
include: os.path.join(RULES, "sickle.rules")
include: os.path.join(RULES, "sra_to_fastq_split.rules")
include: os.path.join(RULES, "subread_align.rules")
include: os.path.join(RULES, "subread_featureCounts.rules")
include: os.path.join(RULES, "subread_index.rules")

# Python includes
PYTHON = os.path.join(GENEREG_LIB, "scripts/python_lib")
include: os.path.join(PYTHON, "util.py")

# Define verbosity
if not ("verbosity" in config.keys()):
    config["verbosity"] = 0
verbosity = int(config["verbosity"])

#================================================================#
#                      Data & wildcards                          #
#================================================================#

# Raw data
if not (("dir" in config.keys()) and ("reads_source" in config["dir"].keys())):
    sys.exit("The parameter config['dir']['reads_source'] should be specified in the config file.")
else:
    READS = config["dir"]["reads_source"]

if not (("metadata" in config.keys()) and ("strands" in config["metadata"].keys())):
    sys.exit("The strands suffixes (parameter config['metadata']['strands']) should be specified in the config file for a paired ends analysis.")

STRANDS = config["metadata"]["strands"].split()

# Samples
SAMPLES = read_table(config["metadata"]["samples"])
SAMPLE_IDS = SAMPLES.iloc[:,0]

## Design ## TODO 
DESIGN = read_table(config["metadata"]["design"])


## Data & results dir

if not (("dir" in config.keys()) and ("reads_source" in config["dir"].keys())):
    sys.exit("The parameter config['dir']['reads_source'] should be specified in the config file.")

READS = config["dir"]["reads_source"]
if not os.path.exists(READS):
    os.makedirs(READS)

if not ("results" in config["dir"].keys()):
    sys.exit("The parameter config['dir']['results'] should be specified in the config file.")
else:
    RESULTS_DIR = config["dir"]["results"]

if not ("samples" in config["dir"].keys()):
    SAMPLE_DIR = config["dir"]["results"]
else:
    SAMPLE_DIR = config["dir"]["samples"]

if not ("reports" in config["dir"].keys()):
    REPORTS_DIR = config["dir"]["results"]
else:
    REPORTS_DIR = config["dir"]["reports"]

if not ("diffexpr" in config["dir"].keys()):
    DEG_DIR = config["dir"]["results"]
else:
    DEG_DIR = config["dir"]["diffexpr"]



#================================================================#
#                         Workflow                               #
#================================================================#

## Data import 
IMPORT = expand(SAMPLE_DIR + "{samples}/{samples}_{strand}.fastq", samples=SAMPLE_IDS, strand=STRANDS)


## Genome & annotations download
GENOME = config["genome"]["version"]
GENOME_DIR = config["dir"]["genome"] + config["genome"]["version"] + "/"

if not os.path.exists(GENOME_DIR):
    os.makedirs(GENOME_DIR)

GENOME_FASTA = expand(GENOME_DIR + "/" + GENOME + ".fa")
GENOME_ANNOTATIONS = expand(GENOME_DIR + "/" + GENOME + ".{ext}", ext=["gff3", "gtf"])

## Graphics & reports
GRAPHICS = expand(REPORTS_DIR + "{graph}.png", graph=["dag", "rulegraph"])

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(SAMPLE_DIR + "{samples}/{samples}_{strand}_fastqc/{samples}_{strand}_fastqc.html", samples=SAMPLE_IDS, strand=STRANDS)
QC = RAW_QC

#----------------------------------------------------------------#
# Trimming
#----------------------------------------------------------------#

if not (("tools" in config.keys()) and ("trimming" in config["tools"].keys())):
    sys.exit("The parameter config['tools']['trimming'] should be specified in the config file. Empty quotes equal to no trimming.")

TRIMMING_SUFFIXES = ["_" + s for s in config["tools"]["trimming"].split()]

TRIMMING = expand(SAMPLE_DIR + "{samples}/{samples}{trimmer}_{strand}", samples=SAMPLE_IDS, trimmer=TRIMMING_SUFFIXES, strand=STRANDS)
TRIM = expand("{trimming}.fastq", trimming=TRIMMING)

TRIM_QC = expand(SAMPLE_DIR + "{samples}/{samples}{trimmer}_{strand}_fastqc/{samples}{trimmer}_{strand}_fastqc.html", samples=SAMPLE_IDS, trimmer=TRIMMING_SUFFIXES, strand=STRANDS)

QC = RAW_QC + TRIM_QC


#----------------------------------------------------------------#
# Alignment
#----------------------------------------------------------------#

if not (("tools" in config.keys()) and ("mapping" in config["tools"].keys())):
    sys.exit("The parameter config['tools']['mapping'] should be specified in the config file.")

MAPPING_TOOLS = config["tools"]["mapping"].split()
MAPPING_SUFFIXES = ["_" + s for s in MAPPING_TOOLS]

ALIGNMENT=expand(SAMPLE_DIR + "{samples}/{samples}{trimmer}{aligner}", samples=SAMPLE_IDS, aligner=MAPPING_SUFFIXES, trimmer=TRIMMING_SUFFIXES)

INDEX = expand(GENOME_DIR + "{tool}/" + GENOME + ".fa", map=MAPPING_TOOLS)

MAPPING = expand("{alignment}.bam", alignment=ALIGNMENT)

SORTED_BAM = expand("{alignment}_sorted_pos.bam", alignment=ALIGNMENT)
BAM_STATS = expand("{alignment}_bam_stats.txt", alignment=ALIGNMENT)


#----------------------------------------------------------------#
# RNA-seq analysis
#----------------------------------------------------------------#

FEATURE_COUNTS = expand("{alignment}_featureCounts.txt", alignment=ALIGNMENT)

INFER_TRANSCRIPTS = expand("{alignment}_cufflinks/transcripts.gtf", alignment=ALIGNMENT)

if not (("tools" in config.keys()) and ("diffexpr" in config["tools"].keys())):
    sys.exit("The parameter config['tools']['diffexpr'] should be specified in the config file.")

DIFFEXPR_TOOLS = config["tools"]["diffexpr"]

SARTOOLS_TARGETFILE = expand(DEG_DIR + "{map}_SARTools_design.txt", map=MAPPING_TOOLS)

DEG = expand(DEG_DIR + "{map}_{deg}_report.html", map=MAPPING_TOOLS, deg=DIFFEXPR_TOOLS)


## ----------------------------------------------------------------
## Visualization & reports
## ----------------------------------------------------------------

## TODO


GENOME_COVERAGE = expand("{alignment}.bedgraph", alignment=ALIGNMENT)
GENOME_COVERAGE_GZ = expand("{alignment}.bedgraph.gz", alignment=ALIGNMENT)

#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses.
	"""
	input: GRAPHICS, QC, GENOME_ANNOTATIONS, FEATURE_COUNTS, INFER_TRANSCRIPTS, SARTOOLS_TARGETFILE, DEG
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

