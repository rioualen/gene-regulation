"""ChIP-seq workflow including: 
	- normalization & QC: fastq merging, trimming, quality control...
	- formatting rules: bed, bam, sam, narrowPeak, fasta...
	- alignment: Bowtie, BWA
	- peak-calling: SWEMBL, MACS, SPP, HOMER
	- downstream analyses: IDR, sequence purge, peak length
"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import time
from subprocess import call

#================================================================#
#                    Variables                                   #
#================================================================#

configfile: "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/config_claire.json"
workdir: config["dir"]["base"] ## does not work?

READS = config["dir"]["reads_source"]
RESULTSDIR = config["dir"]["results"]	

# Rules dir
RULES = "/home/rioualen/Desktop/workspace/fg-chip-seq/scripts/snakefiles/rules/"

#ALN
ALIGNERS = config['aligners'].split()

# Raw data, non auto
CHIP = config["samples"]["chip"].split()
INPUT = config["samples"]["input"].split()
SAMPLES = CHIP + INPUT

# Data merging
for sample in SAMPLES:
	inp = os.listdir(READS+sample)
	inp2 = []
	for i in inp:
		i = READS+sample+"/"+i
		inp2.append(i)
	out = RESULTSDIR + sample + "/" + sample + "_merged.fastq"
	os.system("ls -1 " + ' '.join(inp2) + " | xargs cat > " + out)

#MERGED = expand(RESULTSDIR + "{samples}/{samples}_merged.fastq", samples=SAMPLES)

# Data trimming
SICKLE_TRIMMING = expand(config['dir']['results'] + "{samples}/{samples}_merged_sickle_se_q" + config['sickle']['threshold'] + ".fastq", samples=SAMPLES)

# Quality control
RAW_QC = expand(config["dir"]["results"] + "{samples}/{samples}_merged_fastqc/", samples=SAMPLES)
TRIMMED_QC = expand(config["dir"]["results"] + "{samples}/{samples}_merged_sickle_se_q" + config['sickle']['threshold'] + "_fastqc/", samples=SAMPLES)

# Mapping Bowtie
BOWTIE_MAPPING = expand(config["dir"]["results"] + "{samples}/{samples}_merged_sickle_se_q" + config["sickle"]["threshold"] + "_bowtie.sam", samples=SAMPLES)
BWA_MAPPING = expand(config["dir"]["results"] + "{samples}/{samples}_merged_sickle_se_q" + config["sickle"]["threshold"] + "_bwa.sam", samples=SAMPLES)

# File conversion
SAM_TO_BAM = expand('{sample}/{sample}/{sample}_{aligner}.bam', sample=SAMPLES, aligner=ALIGNERS)
BAM_TO_BED = expand('{sample}/{sample}/{sample}_{aligner}.bed', sample=SAMPLES, aligner=ALIGNERS)

### List all the raw read files, which will be submitted to quality control
#RAWR_FILES, RAWR_DIRS, RAWR_BASENAMES=glob_multi_dir(SAMPLE_IDS, "*_R*_001.fastq.gz", config["dir"]["reads"], ".fastq.gz")

### Quality control
##RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)
#MERGED_RAWR_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"] + "_fastqc/", zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
##MERGED_RAWR_PREFIXES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_fwd"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS) + expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_merged" + config["suffix"]["reads_rev"], zip, reads=PAIRED_BASENAMES, sample_dir=PAIRED_DIRS)
##MERGED_RAWR_QC = expand("{fastq_prefix}_fastqc", fastq_prefix=MERGED_PREFIXES)

### Trimmed reads
## TRIMMED_SUMMARIES = expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_summary.txt", zip, reads=RAWR_BASENAMES_FWD, sample_dir=RAWR_DIRS_FWD)
## TRIMMED_FILES, TRIMMED_DIRS, TRIMMED_BASENAMES=glob_multi_dir(SAMPLE_IDS, "*_R*_001_trimmed_thr" + config["sickle"]["threshold"] + ".fastq.gz", config["dir"]["reads"], ".fastq.gz")
## TRIMMED_QC=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_fastqc/", zip, reads=TRIMMED_BASENAMES, sample_dir=TRIMMED_DIRS)

### Mapped reads
##MAPPED_FILES=expand(config["dir"]["reads"] + "/{sample_dir}/{reads}_trimmed_thr" + config["sickle"]["threshold"] + "_bowtie.sam", zip, reads=RAWR_BASENAMES, sample_dir=RAWR_DIRS)

rule all: 
    """
    Run all the required analyses
    """
    input: 
		expand(RESULTSDIR + "dag.pdf"), \
		SICKLE_TRIMMING, BOWTIE_MAPPING, RAW_QC, TRIMMED_QC, BWA_MAPPING, SAM_TO_BAM, BAM_TO_BED
    params: qsub=config["qsub"]
    shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"



#rule merge_reads:
#    """Merge the raw short read files (SRR) belonging to the same sample (GSM) before allowing the trimming."""
#	input: READS + "GSM570045/SRR063882.fastq", READS + "GSM570045/SRR063883.fastq"
#    	output: RESULTSDIR + "GSM570045/GSM570045_merged.fastq"
#    	shell: "ls -1 {input} | xargs cat > {output};"

#rule merge_reads_2:
#    """Merge the raw short read files (SRR) belonging to the same sample (GSM) before allowing the trimming."""
#	input: READS + "GSM570046/SRR063884.fastq"
#    	output: RESULTSDIR + "GSM570046/GSM570046_merged.fastq"
#    	shell: "ls -1 {input} | xargs cat > {output};"

#================================================================#
#                         Workflow                               #
#================================================================#

#rule all:
#	"""Run the workflow on Short Read Run data.
#	"""
#	input: 
#		expand(RESULTSDIR + "dag.pdf"), \
#		expand(RESULTSDIR + "rule.pdf"), \
#		expand(RESULTSDIR + "{samples}/{samples}_sickle.fastq", samples=SAMPLES)
##		expand(RESULTSDIR + "{samples}/{samples}_{step}_fastqc.html", samples=SAMPLES, step=QUALITY_STEP), \
##		expand(RESULTSDIR + "{samples}/{samples}_{aligner}.sam", samples=SAMPLES, aligner=ALIGNER)
##		expand(RESULTSDIR + "{samples}/{samples}_{trimmer}.fastq", gsm=GSM_LIST, trimmer="sickle"), \
##		expand(RESULTSDIR + "{samples}/{samples}_{trimmer}_fastqc.html", gsm=GSM_LIST, trimmer="sickle")
###		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_{oligo}.txt", chip=CHIP, inp=INPUT, oligo=OLIGO_LENGTH, aligner=ALIGNER, caller=PEAK_CALLER), \			
##		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}_length.png", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER)
##		expand("results/{chip}_vs_{inp}/{chip}_{aligner}_{caller}.bed", chip=CHIP, inp=INPUT, aligner=ALIGNER, caller=PEAK_CALLER)
##		expand("results/IDR/idr.done")

#================================================================#
#                         Includes                               #
#================================================================#

#include: RULES + "bed_to_fasta.rules"
include: RULES + "bowtie.rules"
include: RULES + "bwa.rules"
#include: RULES + "convert_bam_to_bed.rules"
#include: RULES + "convert_sam_to_bam.rules"
#include: RULES + "count_oligo.rules"
include: RULES + "fastqc.rules"
include: RULES + "flowcharts.rules"
#include: RULES + "homer.rules"
#include: RULES + "idr.rules"
#include: RULES + "macs14.rules"
#include: RULES + "merge.rules"
#include: RULES + "narrowpeak_to_bed.rules"
#include: RULES + "peak_length.rules"
#include: RULES + "purge_sequence.rules"
include: RULES + "sickle_se.rules"
#include: RULES + "sorted_bam.rules"
#include: RULES + "spp.rules"
#include: RULES + "swembl.rules"
#include: RULES + "trimming.rules"

