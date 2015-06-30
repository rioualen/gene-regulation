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

GSM_LIST = read_analysis_table_lib.get_gsm_list(config["files"]["analyses"])

include: "rules/flowcharts.rules"
include: "rules/rsync.rules"
include: "rules/fastqc_lucie.rules"
include: "rules/sickle_se.rules"
include: "rules/subread_mapping.rules"


rule all:
    input:expand(config["dir"]["results"] + "{dataset}/{dataset}_fastqc/", dataset = GSM_LIST), \
    expand(config["dir"]["results"] + "{dataset}/{dataset}_trimmed_fastqc/", dataset = GSM_LIST), \
    expand(config["dir"]["results"] + "{dataset}/{dataset}_{aligneur}.bam", dataset = GSM_LIST, aligneur= "subread")