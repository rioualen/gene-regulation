"""Generic work flow for the analysis of ChIP-seq data for the binding
of transcription factors.


This workflow performs the following treatments: 

 - read mapping
 - peak-calling with alternate peak-calling programs
 - motif discovery

The details are specified in a yaml-formatted configuration file.

Usage: 
    snakemake -p  -c "qsub {params.qsub}" -j 12 \
        -s scripts/snakefiles/workflows/factor_workflow.py \
        --configfile path/to/specific/config_file.yml \
        [targets]

Flowcharts:
    snakemake -p -s scripts/snakefiles/workflows/factor_workflow.py \
        --configfile path/to/specific/config_file.yml \
        --force flowcharts

Reference genome:	-
Sequencing type: 	single end

Author: 		Claire Rioualen, Jacques van Helden
Contact: 		claire.rioualen@inserm.fr

"""

#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R
import os
import sys
import datetime
import pandas as pd

## Config
configfile: "examples/GSE20870/GSE20870.yml"

#================================================================#
#                    Check mandatory parameters
#================================================================#

# Define verbosity
if not ("verbosity" in config.keys()):
    config["verbosity"] = 0
verbosity = int(config["verbosity"])

#================================================================#
#                         Includes                               #
#================================================================#
if not ("dir" in config.keys()) & ("fg_lib" in config["dir"].keys()) :
    sys.exit("The parameter config['dir']['fg_lib'] should be specified in the config file.")
FG_LIB = os.path.abspath(config["dir"]["fg_lib"])
RULES = os.path.join(FG_LIB, "scripts/snakefiles/rules")
PYTHON = os.path.join(FG_LIB, "scripts/python_lib")

include: os.path.join(PYTHON, "util.py")

include: os.path.join(RULES, "annotation_download.rules")
include: os.path.join(RULES, "bam_by_name.rules")
include: os.path.join(RULES, "bam_by_pos.rules")
include: os.path.join(RULES, "bam_to_bed.rules")
include: os.path.join(RULES, "bam_stats.rules")
include: os.path.join(RULES, "bowtie_index.rules")
include: os.path.join(RULES, "bowtie_se.rules")
include: os.path.join(RULES, "bowtie2_index.rules")
include: os.path.join(RULES, "bowtie2_se.rules")
include: os.path.join(RULES, "bPeaks.rules")
include: os.path.join(RULES, "bwa_index.rules")
include: os.path.join(RULES, "bwa_se.rules")
include: os.path.join(RULES, "count_reads.rules")
include: os.path.join(RULES, "fastqc.rules")
include: os.path.join(RULES, "flowcharts.rules")
include: os.path.join(RULES, "genome_coverage_bedgraph.rules")
include: os.path.join(RULES, "genome_download.rules")
include: os.path.join(RULES, "getfasta.rules")
include: os.path.join(RULES, "get_chrom_sizes.rules")
include: os.path.join(RULES, "gzip.rules")
include: os.path.join(RULES, "homer.rules")
#include: os.path.join(RULES, "import_fastq.rules")
include: os.path.join(RULES, "igv_session.rules")
include: os.path.join(RULES, "macs2.rules")
include: os.path.join(RULES, "macs14.rules")
include: os.path.join(RULES, "peak_motifs.rules")
include: os.path.join(RULES, "sickle_se.rules")
include: os.path.join(RULES, "spp.rules")
include: os.path.join(RULES, "swembl.rules")
include: os.path.join(RULES, "sam_to_bam.rules")
include: os.path.join(RULES, "sra_to_fastq.rules")

ruleorder: bam_by_pos > sam_to_bam
ruleorder: bam_by_name > sam_to_bam
#ruleorder: import_fastq > sickle_se 
#================================================================#
#                      Data & wildcards                          #
#================================================================#

# Raw data
READS = config["dir"]["reads_source"]

# Samples
SAMPLES = read_table(config["files"]["samples"])
SAMPLE_IDS = SAMPLES.iloc[:,0]

## Design
DESIGN = read_table(config["files"]["design"])
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


#================================================================#
#                         Workflow                               #
#================================================================#

## Data import 

IMPORT = expand(SAMPLE_DIR + "{samples}/{samples}.fastq", samples=SAMPLE_IDS)

# Genome
GENOME = config["genome"]["version"]
GENOME_DIR = config["dir"]["genomes"] + config["genome"]["version"]
GENOME_DIR = config["dir"]["genomes"] + config["genome"]["version"]

GENOME_FASTA = expand(config["dir"]["genomes"] + config["genome"]["version"] + "/" + config["genome"]["version"] + ".fa")
GENOME_ANNOTATIONS = expand(config["dir"]["genomes"] + config["genome"]["version"] + "/" + config["genome"]["version"] + ".gff3")

### Graphics & reports
GRAPHICS = expand(RESULTS_DIR + "dag.pdf")
#REPORT = expand(RESULTS_DIR + "report.html")

#----------------------------------------------------------------#
# Quality control
#----------------------------------------------------------------#

RAW_QC = expand(SAMPLE_DIR + "{samples}/{samples}_fastqc/{samples}_fastqc.html", samples=SAMPLE_IDS)


#----------------------------------------------------------------#
# Trimming
#----------------------------------------------------------------#

TRIMMER="sickle-se-q" + config["sickle"]["threshold"]
TRIMMING=expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}", samples=SAMPLE_IDS, trimmer=TRIMMER)
TRIM = expand(SAMPLE_DIR + "{trimming}.fastq", trimming=TRIMMING)

TRIM_QC = expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}_fastqc/{samples}_{trimmer}_fastqc.html", samples=SAMPLE_IDS, trimmer=TRIMMER)
QC = RAW_QC + TRIM_QC


#----------------------------------------------------------------#
# Alignment
#----------------------------------------------------------------#


ALIGNER=["bwa"]
ALIGNMENT=expand(SAMPLE_DIR + "{samples}/{samples}_{trimmer}_{aligner}", samples=SAMPLE_IDS, aligner=ALIGNER, trimmer=TRIMMER)

INDEX = expand(config["dir"]["genomes"] + config["genome"]["version"] + "/{aligner}/" + config["genome"]["version"] + ".fa", aligner=ALIGNER)

MAPPING = expand("{alignment}.sam", alignment=ALIGNMENT)

BAM_STATS = expand("{alignment}_bam_stats.txt", alignment=ALIGNMENT)

GENOME_COVERAGE = expand("{alignment}.bedgraph", alignment=ALIGNMENT)
GENOME_COVERAGE_GZ = expand("{alignment}.bedgraph.gz", alignment=ALIGNMENT)



# Sort mapped reads

## Why not work ?
SORTED_BY_POS = expand(SAMPLE_DIR + "{alignment}_sorted_pos.bam", alignment=ALIGNMENT)
SORTED_BY_NAME = expand(SAMPLE_DIR + "{alignment}_sorted_name.bam", alignment=ALIGNMENT)
#BAM_READNB = expand(RESULTS_DIR + "{alignment}_sorted_pos_bam_readnb.txt", alignment=ALIGNMENT)
SORTED_READS_BED = expand(SAMPLE_DIR + "{alignment}_sorted_pos.bed", alignment=ALIGNMENT)
#BED_FEAT_COUNT = expand(RESULTS_DIR + "{alignment}_sorted_pos_bed_nb.txt", alignment=ALIGNMENT)


#TDF = expand(RESULTS_DIR + "{alignment}_sorted_pos.tdf", alignment=ALIGNMENT)

# ----------------------------------------------------------------
# Peak-calling
# ----------------------------------------------------------------


##TODO check if params are defined; if not, set them.
PEAKCALLER=[
#    "homer-fdr" + config["homer"]["fdr"] + "_peaks", 
#    "macs2-qval" + config["macs2"]["qval"], 
#    "swembl-R" + config["swembl"]["R"],
    "macs14-pval" + config["macs14"]["pval"],
    "spp-fdr" + config["spp"]["fdr"],
    "bPeaks-log" + config["bPeaks"]["log2FC"]
]

PEAKCALLING=expand(expand(PEAKS_DIR + "{treat}_vs_{control}/{{peakcaller}}/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER, trimmer=TRIMMER)

PEAKS = expand("{peakcalling}.bed", peakcalling=PEAKCALLING)

# ----------------------------------------------------------------
# Peak analysis
# ----------------------------------------------------------------


MOTIFS=expand(expand("{treat}_vs_{control}/{{peakcaller}}/peak-motifs/{treat}_vs_{control}_{{trimmer}}_{{aligner}}_{{peakcaller}}_peak-motifs_synthesis", zip, treat=TREATMENT, control=CONTROL), peakcaller=PEAKCALLER, aligner=ALIGNER, trimmer=TRIMMER)


GET_FASTA = expand(PEAKS_DIR + "{peakcalling}.fasta", peakcalling=PEAKCALLING)
PEAK_MOTIFS = expand(PEAKS_DIR + "{motifs}.html", motifs=MOTIFS)

# ----------------------------------------------------------------
# Visualization
# ----------------------------------------------------------------

VISU = expand(PEAKS_DIR + "igv_session.xml")

#================================================================#
#                        Rule all                                #
#================================================================#

rule all: 
	"""
	Run all the required analyses.
	"""
	input: GRAPHICS, BAM_STATS, PEAKS, QC, GENOME_COVERAGE_GZ, GENOME_ANNOTATIONS, VISU#PEAK_MOTIFS#, CHROM_SIZES, PEAKS, TDFRAW_QC, MAPPING, PEAKS, IMPORT, INDEX, PEAKS, 
	params: qsub=config["qsub"]
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"

#================================================================#
#                        IGV stuff                               #
#================================================================#

#filename = PEAKS_DIR + "IGV_session.xml"

#file = open(filename, "w")



#file.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
#file.write('<Session genome="' + GENOME_FASTA[0] + '" hasGeneTrack="true" hasSequenceTrack="true" locus="all" path="' + filename + '" version="8">\n')

#file.write('    <Resources>\n')
#for i in PEAKS:
#    file.write('        <Resource path="' + i + '"/>\n')
#for i in GENOME_COVERAGE_GZ:
#    file.write('        <Resource path="' + i + '"/>\n')
#file.write('        <Resource path="' + GENOME_ANNOTATIONS[0] + '"/>\n')
#file.write('    </Resources>\n\n')

#file.write('    <Panel height="519" name="DataPanel" width="1901">\n')
#for i in PEAKS:
#    file.write('        <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,153,255" colorScale="ContinuousColorScale;0.0;52.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" id="' + i + '" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">\n')
#    file.write('                <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="52.0" minimum="0.0" type="LINEAR"/>\n')
#    file.write('        </Track>\n')
#file.write('    </Panel>\n\n')


#file.write('    <Panel height="259" name="AlignmentPanel" width="1901">\n')
#for i in GENOME_COVERAGE_GZ:
#    file.write('        <Track height="50" clazz="org.broad.igv.track.DataSourceTrack" color="113,35,30" displayMode="COLLAPSED" featureVisibilityWindow="-1" id="' + i + '" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="max">\n')
#    file.write('            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="200.0" minimum="0.0" type="LINEAR"/>\n')
#    file.write('        </Track>\n')
#file.write('    </Panel>\n\n')

#file.write('    <Panel height="120" name="FeaturePanel" width="1901">\n')
#file.write('        <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;235.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="/data/genomes/sacCer3/sacCer3.gff3" name="sacCer3.gff3" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">\n')
#file.write('            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="235.0" minimum="0.0" type="LINEAR"/>\n')
#file.write('        </Track>\n')
#file.write('    </Panel>\n\n')


#file.write('    <PanelLayout dividerFractions="0.332484076433121"/>\n')
#file.write('    <HiddenAttributes>\n')
#file.write('        <Attribute name="NAME"/>\n')
#file.write('        <Attribute name="DATA FILE"/>\n')
#file.write('        <Attribute name="DATA TYPE"/>\n')
#file.write('    </HiddenAttributes>\n')
#file.write('</Session>\n')




#file.close()
