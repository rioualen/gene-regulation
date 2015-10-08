"""Looking for parB sites in P.aeruginosa.

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/Scerevisiae-Pho4.py \
        --force flowcharts

Organism: 		Pseudomonas aeruginosa
Reference genome:	-
Sequencing type: 	single end
Data source: 		Gene Expression Omnibus

Author: 		Claire Rioualen, Jacques van Helden
Contact: 		claire.rioualen@inserm.fr
"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time#rm?
import datetime
import pandas as pd



## Config
configfile: "scripts/snakefiles/workflows/Paeruginosa_local.json"
workdir: config["dir"]["base"]
verbosity = int(config["verbosity"])

#================================================================#
#                         Includes                               #
#================================================================#

FG_LIB = os.path.abspath(config["dir"]["fg_lib"])
RULES = os.path.join(FG_LIB, "scripts/snakefiles/rules")
PYTHON = os.path.join(FG_LIB, "scripts/snakefiles/python_lib")

include: os.path.join(PYTHON, "util.py")
include: os.path.join(RULES, "util.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "bwa_index.rules")
include: os.path.join(RULES, "bwa_se.rules")
include: os.path.join(RULES, "bed_to_fasta.rules")
include: os.path.join(RULES, "convert_bam_to_bed.rules")
include: os.path.join(RULES, "count_oligo.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "flowcharts.rules")
include: os.path.join(RULES, "homer.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "peak_length.rules")
#include: os.path.join(RULES, "peak_motifs.rules")
#include: os.path.join(RULES, "purge_sequence.rules")
include: os.path.join(RULES, "sickle_se.rules")
#include: os.path.join(RULES, "sorted_bam.rules")
include: os.path.join(RULES, "spp.rules")
#include: os.path.join(RULES, "sra_to_fastq.rules")
include: os.path.join(RULES, "swembl.rules")

#================================================================#
#                      Data & config                             #
#================================================================#

# Raw data
READS = config["dir"]["reads_source"]

# Samples
SAMPLES = read_table(config["files"]["samples"], verbosity=verbosity)
SAMPLE_IDS = SAMPLES.iloc[:,0] ## First column MUST contain the sample ID


## Design
DESIGN = read_table(config["files"]["design"], verbosity=verbosity)
TREATMENT = DESIGN.iloc[:,0]
CONTROL = DESIGN.iloc[:,1]

## Ref genome
GENOME = config["genome"]["genome_version"]

## Results dir
RESULTS_DIR = config["dir"]["results"]
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

## Programs
TRIMMING = "sickle-se-q" + config["sickle"]["threshold"]
ALIGNER = "bwa"

#================================================================#
#                         Workflow                               #
#================================================================#

## Data import & merging.

## À revoir, rsync output du texte -> plante flowcharts...
#os.system("rsync -rupltv --exclude '*.tab' " + READS + " " + RESULTS_DIR)

## Graphics & reports
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")
REPORT = expand(RESULTS_DIR + "report.html")

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(RESULTS_DIR + "{samples}/{samples}_fastqc/", samples=SAMPLE_IDS)

RAW_READNB = expand(RESULTS_DIR + "{samples}/{samples}_fastq_readnb.txt", samples=SAMPLE_IDS)

#----------------------------------------------------------------#
# Trimming (if need be)
#----------------------------------------------------------------#

#TRIMMED_READS_SICKLE = expand(RESULTS_DIR + "{samples}/{samples}_sickle-se-q" + config["sickle"]["threshold"] + ".fastq", samples=SAMPLE_IDS)
#TRIMMED_QC = expand(RESULTS_DIR + "{samples}/{samples}_sickle-se-q" + config["sickle"]["threshold"] + "_fastqc/", samples=SAMPLE_IDS)

#----------------------------------------------------------------#
# Mapping
#----------------------------------------------------------------#

BWA_INDEX = expand(config["dir"]["genome"] + "{genome}/BWAIndex/{genome}.fa.bwt", genome=GENOME)
BWA_MAPPING = expand(RESULTS_DIR + "{samples}/{samples}_{aligner}.sam", samples=SAMPLE_IDS, aligner=ALIGNER)

# Sorted and converted reads (bam, bed)
SORTED_MAPPED_READS_BWA = expand(RESULTS_DIR + "{sample}/{sample}_{aligner}_sorted_pos.bam", sample=SAMPLE_IDS, aligner=ALIGNER)
BAM_READNB = expand(RESULTS_DIR + "{samples}/{samples}_{aligner}_sorted_pos_bam_readnb.txt", samples=SAMPLE_IDS, aligner=ALIGNER)
SORTED_READS_BED = expand(RESULTS_DIR + "{sample}/{sample}_{aligner}_sorted_pos.bed", sample=SAMPLE_IDS, aligner=ALIGNER)
BED_READNB = expand(RESULTS_DIR + "{samples}/{samples}_{aligner}_sorted_pos_bed_nb.txt", samples=SAMPLE_IDS, aligner=ALIGNER)

# ----------------------------------------------------------------
# Peak-calling
# ----------------------------------------------------------------

PEAKS_MACS2 = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}_{{aligner}}_macs2_peaks.narrowPeak",
               zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)

PEAKS_SWEMBL = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/swembl/{treat}_vs_{ctrl}_{{aligner}}_swembl-R" + config["swembl"]["R"] + ".bed", 
               zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)

PEAKS_SPP = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/spp/{treat}_vs_{ctrl}_{{aligner}}_spp-fdr" + config["spp"]["fdr"] + ".narrowPeak", 
               zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)

PEAKS_HOMER = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}_{{aligner}}_homer_peaks.bed", 
              zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)

PEAKS = PEAKS_MACS2 + PEAKS_SWEMBL + PEAKS_HOMER + PEAKS_SPP

# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------

# File conversion / fetching fasta
FETCH_MACS2_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}_{{aligner}}_macs2_peaks.fasta", zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)
FETCH_SWEMBL_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/swembl/{treat}_vs_{ctrl}_{{aligner}}_swembl-R0.01.fasta", zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)
FETCH_SPP_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/spp/{treat}_vs_{ctrl}_{{aligner}}_spp.fasta", zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)
FETCH_HOMER_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}_{{aligner}}_homer_peaks.fasta", zip, treat=TREATMENT, ctrl=CONTROL), aligner=ALIGNER)

FETCH_PEAKS = FETCH_MACS2_PEAKS + FETCH_HOMER_PEAKS

## Sequence purge
#PURGE_MACS2_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)
#PURGE_SWEMBL_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/swembl/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_swembl-R0.01_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)
#PURGE_SPP_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/spp/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_spp_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)
#PURGE_HOMER_PEAKS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)

#PURGE_PEAKS = PURGE_MACS2_PEAKS + PURGE_HOMER_PEAKS

### Working til here ##

## Oligo analysis ## BUG oligo.sam...
#OLIGO = config['oligo_analysis']['count_oligo'].split()
#MACS2_OLIGO = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged_oligo{{oligo}}.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER, oligo=OLIGO)
#HOMER_OLIGO = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged_oligo{{oligo}}.fasta", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER, oligo=OLIGO)

## Peaks length ## TODO when purge works -> vmatch issue
#MACS2_PEAKS_LENGTH = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged_length.png", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)
#HOMER_PEAKS_LENGTH = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged_length.png", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)

## Peak-motifs
#MACS2_MOTIFS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/macs2/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_macs2_peaks_purged.html", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)
#HOMER_MOTIFS = expand(expand(RESULTS_DIR + "{treat}_vs_{ctrl}/homer/{treat}_vs_{ctrl}{{trimming}}_{{aligner}}_homer_peaks_purged.html", zip, treat=TREATMENT, ctrl=CONTROL), trimming=TRIMMING, aligner=ALIGNER)


## Rule Order: pipeline still works when commented?...
#ruleorder: macs2 > bam_to_bed > sam2bam
#ruleorder: count_reads_bam > sam2bam
#ruleorder: bed_to_fasta > purge_sequence


#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses
	"""
##	input: GRAPHICS, IMPORT, TRIMMED_READS_SICKLE, TRIMMED_QC, RAW_QC, MAPPED_READS_BWA, RAW_READNB, BAM_READNB, BED_READNB, PEAKS_MACS2, FETCH_MACS2_PEAKS, PURGE_MACS2_PEAKS #redundant for flowcharts
	input: GRAPHICS, RAW_QC, PEAKS
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"






#================================================================#
#                          Report                                #
#================================================================#


##----------------------------------------------------------------#
## Build the report (including DAG and rulegraph flowcharts).
#from snakemake.utils import report

## Bulleted list of samples for the report
#SAMPLE_IDS_OL=report_numbered_list(SAMPLE_IDS)
#RAW_READS_OL=report_numbered_list(IMPORT)
#TRIMMED_READS_OL=report_numbered_list(TRIMMED_READS_SICKLE)
#RAW_QC_OL=report_numbered_list(RAW_QC)
#TRIMMED_QC_OL=report_numbered_list(TRIMMED_QC)

#MAPPED_SAM_OL=report_numbered_list(MAPPED_READS_BWA)
#MAPPED_BAM_SORTED_OL=report_numbered_list(SORTED_MAPPED_READS_BWA)
#MAPPED_BED_SORTED_OL=report_numbered_list(SORTED_READS_BED)

#RAW_QC_OL=report_numbered_list(RAW_QC)
#TRIMMED_QC_OL=report_numbered_list(TRIMMED_QC)

#PEAKFILES_MACS2_OL=report_numbered_list(PEAKS_MACS2)

#	input: GRAPHICS, IMPORT, TRIMMED_READS_SICKLE, TRIMMED_QC, RAW_QC, MAPPED_READS_BWA, RAW_READNB, BAM_READNB, BED_READNB, PEAKS_MACS2, FETCH_MACS2_PEAKS, PURGE_MACS2_PEAKS #redundant for flowcharts

NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

#rule report:
#    """
#    Generate a report with the list of datasets + summary of the results.
#    """
#    input:  dag=config["dir"]["reports"] + "dag.pdf", \
#            dag_png=config["dir"]["reports"] + "dag.png", \
#            rulegraph=config["dir"]["reports"] + "rule.pdf", \
#            rulegraph_png=config["dir"]["reports"] + "rule.png"
#    output: html=config["dir"]["reports"] + "report.html"
#    run:
#        report("""
#        ===========================================
#        ChIP-seq analysis - S.cerevisiae Pho4 study
#        ===========================================
#        
#        :Date:                 {NOW}
#        :Project:              S cerevisiae
#        :Analysis workflow:    Claire Rioualen
#        
#        Contents
#        ========
#        
#        - `Flowcharts`_
#        - `Datasets`_
#             - `Samples`_
#             - `Raw reads`_
#             - `Trimmed`_
#             - `Mapped`_
#             - `Peaks`_
#             - `QC reports`_

#        -----------------------------------------------------

#        Flowcharts
#        ==========

#        - Sample treatment: dag_
#        - Workflow: rulegraph_

#        .. image:: rulegraph.png

#        -----------------------------------------------------

#        Datasets
#        ========
#        
#        Samples
#        -------

#        {SAMPLE_IDS_OL} 

#        Raw reads 
#        ---------

#        {RAW_READS_OL}

#        Trimmed
#        -------

#        {TRIMMED_READS_OL}

#        Mapped
#        ------

#        Sam format (uncompressed)

#        {MAPPED_SAM_OL}

#        Bam format (sorted by positions)

#        {MAPPED_BAM_SORTED_OL}

#        Bed format (sorted by positions)

#        {MAPPED_BED_SORTED_OL}

#        Peaks
#        -----

#        Macs2 peaks

#        {PEAKFILES_MACS2_OL}

#        QC reports
#        ----------

#        {RAW_QC_OL}

#        {TRIMMED_QC_OL}

#        -----------------------------------------------------

#        """, output.html, metadata="Claire Rioualen (claire.rioualen@inserm.fr)", **input)




