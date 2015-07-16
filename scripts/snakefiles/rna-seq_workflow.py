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
configfile: "scripts/snakefiles/config_rna-seq.json"
WDIR = config["dir"]["base"]
workdir: WDIR


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
include: "rules/fastqc.rules"
include: "rules/sickle_se.rules"
include: "rules/subread_mapping.rules"
include: "rules/featurecounts.rules"

COUNT_FILES = expand(config["dir"]["results"] + "{dataset}/{dataset}_{trimmed}_{aligneur}_featurecounts.tab", dataset = GSM_LIST, aligneur= "subread", trimmed = "sickle_se_q" + config["sickle"]["threshold"] )
PARAMS_R = config["dir"]["results"] + "DEG/sickle_se_q" + config["sickle"]["threshold"] + "_subread_featurecounts_params.R"
ALL_COUNTS = config["dir"]["results"] + "DEG/sickle_se_q" + config["sickle"]["threshold"] + "_subread_featurecounts_allcounts.tab"

include: "rules/HTseq_allcount_params.rules"

RESULTS_EDGER = expand(config["dir"]["results"] + "DEG/{cond_1}_vs_{cond_2}/{cond_1}_VS_{cond_2}_{trimmed}_{aligneur}_{count}" + config["edgeR"]["software"] +".tab", zip, cond_1 = config["Diff_Exp"]["cond1"], cond_2 = config["Diff_Exp"]["cond2"], aligneur= "subread", trimmed = "sickle_se_q" + config["sickle"]["threshold"], count = "featurecounts")

include: "rules/edgeR.rules"

RESULTS_DESEQ2 = expand(config["dir"]["results"] + "DEG/{cond_1}_vs_{cond_2}/{cond_1}_VS_{cond_2}_{trimmed}_{aligneur}_{count}" + config["DESeq2"]["software"] +".tab", zip, cond_1 = config["Diff_Exp"]["cond1"], cond_2 = config["Diff_Exp"]["cond2"], aligneur= "subread", trimmed = "sickle_se_q" + config["sickle"]["threshold"], count = "featurecounts")
include: "rules/DESeq2.rules"


rule all:
    input:expand(config["dir"]["results"] + "{dataset}/{dataset}_fastqc/", dataset = GSM_LIST), \
    expand(config["dir"]["results"] + "{dataset}/{dataset}_sickle_se_q" + config["sickle"]["threshold"] + "_fastqc/", dataset = GSM_LIST), \
    expand(config["dir"]["results"] + "{dataset}/{dataset}_{trimmed}_{aligneur}_featurecounts.tab", dataset = GSM_LIST, aligneur= "subread", trimmed = "sickle_se_q" + config["sickle"]["threshold"] ), \
    RESULTS_EDGER, RESULTS_DESEQ2