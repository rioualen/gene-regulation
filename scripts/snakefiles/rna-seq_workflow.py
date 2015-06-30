"""Generic workflow for RNA-seq analysis.

This workflow combines a series of rules generally required to ensure
the conventional treatment of RNA-seq transcriptomic data. 

- read mapping
- quality control
- estimation of gene abundance
- detection of differentially expressed genes (DEG)
- functional enrichment of DEG
- ...
"""

WDIR = "/home/lkhamvongsa/dr-chip-rna-seq"
workdir: WDIR
configfile: "scripts/snakefiles/config_rna-seq.json"

from snakemake.utils import R
import os
import sys
path = "scripts/python-scripts/"
sys.path.append(path) # For find my module "read_analysis_table_lib"
import read_analysis_table_lib # read_analysis_table_lib is a module where you can find two functions who give you several lists for peak-callers softwares

## Usage: snakemake -c "qsub {params.qsub}" -j 10

"""
    ===========================
    Definition of several lists
    ===========================

"""
ANALYSIS_TABLE = "data/rna-seq/analysis_description.tab"
GSM_LIST = read_analysis_table_lib.get_gsm_list(ANALYSIS_TABLE)

include: "rules/fastqc_lucie.rules"
include: "rules/trimming.rules"
include: "rules/subread_mapping.rules"

rule all:
    input:expand( config["results_directory"] + "/{dataset}/{dataset}_fastqc/", dataset = GSM_LIST), \
    expand(config["results_directory"] + "/{dataset}/{dataset}_{aligneur}.bam", dataset = GSM_LIST, aligneur= "subread")