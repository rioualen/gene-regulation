"""Workflow for the analysis of transcriptomic RNA-seq data in
Desulfovibrio desulfuricans.

Reference genome: Desulfovibrio desulfuricans ATCC_27774 uid59213

Collaboration with Alain Dolla
Data from : Jeff Cole

"""

#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
configfile: "scripts/snakefiles/workflows/no-stress_analysis_config.json"


#================================================================#
#     Global variables                                           #
#================================================================#

# Usage: snakemake -c "qsub {params.qsub}" -j 25
# workdir:"/home/desulfo-no" # Root directory for the project. Should be adapted for porting.
# workdir: "/home/desulfo-no/no-stress_project/data/1258-BRM" # Root directory for the project. Should be adapted for porting.
# workdir:"/home/no-stress/Documents/no-stress_project/data/1258-BRM" # Local Root directory for the project. Should be adapted for porting.
workdir: os.getcwd() # Local Root directory for the project. Should be adapted for porting.
DATASETS = "N1_1 N1_2 N2_1 N2_2 N4_1 N4_2 S1_1 S1_2 S4_1 S4_2 S5_1 S5_2 NN2_1 NN2_2 NN4_1 NN4_2 NN5_1 NN5_2 SN1_1 SN1_2 SN2_1 SN2_2 SN5_1 SN5_2".split() # list of files
DATA_DIRS = "N1 N1 N2 N2 N4 N4 S1 S1 S4 S4 S5 S5 NN2 NN2 NN4 NN4 NN5 NN5 SN1 SN1 SN2 SN2 SN5 SN5".split() # list all directories

# DATASETS = "N1_1 N1_2 NN2_1 NN2_2".split() # list of files
# DATA_DIRS = "N1 N1 NN2 NN2".split() # list the directories for each file



QSUB_PARAM = config["qsub_parameters"]

# Variable for fastqc
OPTIONS_FASTQC = config["fastqc_parameters"]["other_options"]

# Variables for sickle (trimming)
THRESHOLD = config["sickle_parameters"]["threshold"]
SEQ_TYPE = config["sickle_parameters"]["seq_type"]
OPTIONS_SICKLE = config["sickle_parameters"]["other_options"]

# Variables for bowtie2 (alignement)
INDEX = config["bowtie2_parameters"]["bowtie_index"]
MAX_MISMATCHES = config["bowtie2_parameters"]["max_mismatches"]
OPTIONS_BOWTIE2 = config["bowtie2_parameters"]["other_options"]

# Variables for HTseq (count table)
HT_TYPE = config["htseq_parameters"]["ht_type"]
STRANDED = config["htseq_parameters"]["stranded"]
MINAQUAL = config["htseq_parameters"]["minaqual"]
IDATTR = config["htseq_parameters"]["idattr"]
GFF_FILE = config["htseq_parameters"]["gff_file"]
OPTIONS_COUNT = config["htseq_parameters"]["other_options"]


#================================================================#
#                         Includes                               #
#================================================================#

include: "scripts/snakefiles/rules/flowcharts.rules"
include: "scripts/snakefiles/rules/fastqc.rules"
include: "scripts/snakefiles/rules/trimming.rules"
include: "scripts/snakefiles/rules/bowtie2.rules"
include: "scripts/snakefiles/rules/convert_sam_to_bam.rules"
include: "scripts/snakefiles/rules/sorted_bam.rules"
include: "scripts/snakefiles/rules/count_table.rules"
# include: "scripts/snakefiles/rules/index_bam.rules" # to devlopp if needed for IGV

#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    """Run workflow for each replica of each experience"""
    input:
        expand(config["data_root_dir"] + "{data_dir}/{dataset}_fastqc/", zip, dataset=DATASETS, data_dir=DATA_DIRS), \
        expand(config["data_root_dir"] + "{data_dir}/{dataset}_trimmed_thr20_fastqc/", zip, dataset=DATASETS, data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_trimmed_thr" + THRESHOLD + ".fastq.gz", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + MAX_MISMATCHES + ".sam", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + MAX_MISMATCHES + ".bam", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + MAX_MISMATCHES + "_sn.bam", data_dir=DATA_DIRS), \
        expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + MAX_MISMATCHES + "_count.txt", data_dir=DATA_DIRS), \
        # expand(config["data_root_dir"] + "{data_dir}/{data_dir}_bowtie2_mm" + MAX_MISMATCHES + "_index.bam", data_dir=DATA_DIRS), \ to devlopp if needed for IGV
        
    shell: "echo 'job done'"
