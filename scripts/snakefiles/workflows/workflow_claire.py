#================================================================#
#                        Imports                                 #
#================================================================#

from snakemake.utils import R

"""Workflow applying various motif discovery algorithms in a peak set from ChIP-seq.
"""

#workdir: "~/Desktop/workspace/project2/"	## does not work yet

#================================================================#
#     Global variables (to include in config file instead?)      #
#================================================================#

#QSUB_PARAM = " -V -m a "
GENOME="mm9"
OLIGO_LENGTH = "1 2 4 6".split()
FACTORS="CEBPA HNF4A".split()
SWEMBL_R="0.05 0.02 0.01".split()

#PEAKSDIR = "bedfiles/"
#PEAKS, = glob_wildcards(PEAKSDIR + "{peaks}.bed")

RESULTSDIR = "results/"

#================================================================#
#                         Includes                               #
#================================================================#

include: "fastqc.rules"	
include: "bed_to_fasta.rules"	
include: "purge_sequence.rules"
include: "count_oligo.rules"
include: "peak_length.rules"
include: "flowcharts.rules"


#================================================================#
#                         Workflow                               #
#================================================================#

rule all:
	"""Run the workflow on each peak set.
	"""
 	## what if no output file? 
	## keep all expand even if rules outputs are chained ? -> ! check pb with flowcharts redundancy
	## input: expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_purge.fasta", factor=FACTORS, swembl_r=SWEMBL_R) 
	input: 
		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof.fasta", factor=FACTORS, swembl_r=SWEMBL_R), \
		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_purge.fasta", factor=FACTORS, swembl_r=SWEMBL_R), \
		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_{oligo}.txt", factor=FACTORS, swembl_r=SWEMBL_R, oligo=OLIGO_LENGTH), \
		expand("results/SWEMBL_mmus_{factor}_vs_mmus_Input_peaks_R{swembl_r}_nof_length.png", factor=FACTORS, swembl_r=SWEMBL_R), \
		expand("dag.pdf"), \
		expand("rulegraph.pdf"), \
		expand("{qual}/{reads}_fastqc.html", qual=QUALITY, reads=FASTQ)



	


