# Author Claire Rioualen
# This is a draft


## TODO 
#- chose output location and adapt paths to workdir //!\\
#- check params
#    - mandatory: filename, genome
#    - optional: all other fields (peaks, coverage files, annotation files, etc)
#- check params format: have to be strings or lists of strings, depending
#- fileformats have to be checked:
#    - list of peaks have to be bed-formatted
#    - coverage files can be in bedgraph.gz or tdf (or other?)
#    - genome has to be a fasta/fa file
#    - etc
#- add fields:
#    - transcripts gtf from cufflinks
#    - bam files (! can be very slow to load, relies on bam.bai files)
#    - additional annotation files (like RegulonDB sites or other files mentionned in the config.yml)
#- by default the session will display *all* peakfiles, *all* coverage files, etc. 
#  Later on I need to think of a clean way for the user to customize what to display (only peaks from homer, only bam files from bowtie, whatever)


import pandas as pd
import numpy as np
import os

def igv_session(filename, genome, gtf="", gff3="", peaks="", coverage="", transcripts_gtf=""):

    genome_version = config["genome"]["version"]
    cwd = os.getcwd()

    file = open(filename, "w+")

    file.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    file.write('<Session genome="' + cwd + "/" + genome + '" hasGeneTrack="true" hasSequenceTrack="true" locus="" path="' + filename + '" version="8">\n')

    ## Resource files
    file.write('<Resources>\n')
    if peaks:
        for i in peaks:
            file.write('  <Resource path="../../' + i + '"/>\n')
    if coverage:
        for i in coverage:
            file.write('  <Resource path="../../' + i + '"/>\n')
    if gtf:
        file.write('  <Resource path="../../' + gtf + '"/>\n')
    if gff3:
        file.write('  <Resource path="../../' + gff3 + '"/>\n')
    file.write('</Resources>\n')

    ## Genome annotation panel
    file.write('<Panel name="GeneAnnotPanel" height="60">\n')
    if gff3:
        file.write('  <Track id="../../' + gff3 + '" name="gff3 ' + genome_version + '"  color="153,0,51" fontSize="12" >\n')
        file.write('    <DataRange/>\n')
        file.write('  </Track>\n')
    if gtf:
        file.write('  <Track id="../../' + gtf + '" name="gtf ' + genome_version + '"  color="153,0,51" fontSize="12" >\n')
        file.write('    <DataRange/>\n')
        file.write('  </Track>\n')
    file.write('</Panel>\n')

    ## Genome coverage panel
    if coverage:
        file.write('<Panel name="GenomeCovPanel" height="' + str(len(coverage)*60) + '">\n')
        for bdg in coverage:
            name = bdg.split(sep="/")[-1]
            name = name.split(sep=".")[0]

#            tab = pd.read_table(bdg)
#            cov = tab.iloc[:,3]
#            max = int(np.percentile(cov, 99))

            file.write('  <Track height="50" autoscale="true" id="../../' + bdg + '" name="' + name + '" color="51,153,0" fontSize="12">\n')
            file.write('    <DataRange minimum="0" maximum="200"/>\n')
            file.write('  </Track>\n')
        file.write('</Panel>\n')

    ## Peaks panel
    if peaks:
        file.write('<Panel name="FeaturePanel" height="' + str(len(peaks)*60) + '">\n')
        for bedfile in peaks:
            pc_name = bedfile.split(sep="_")[-1]
            pc_name = os.path.splitext(pc_name)[0]

            sam_name = bdg.split(sep="/")[-1]
            sam_name = sam_name.split(sep=".")[0]

            file.write('  <Track id="../../' + bedfile + '" name="' + pc_name + " " + sam_name + '" color="0,0,178" fontSize="12">\n')
            file.write(    '<DataRange/>\n')
            file.write('  </Track>\n')
        file.write('</Panel>\n')

    file.write('<PanelLayout/>\n')
    file.write('  <HiddenAttributes>\n')
    file.write('    <Attribute name="NAME"/>\n')
    file.write('    <Attribute name="DATA FILE"/>\n')
    file.write('    <Attribute name="DATA TYPE"/>\n')
    file.write('  </HiddenAttributes>\n')
    file.write('</Session>\n')

    file.close()
