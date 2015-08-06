################################################
## Taking as input a result of differential gene 
## expression analysis (e.g. produced by edgeR or DESeq2), 
## this script runs functional enrichment analysis with the 
## BioConductor package clusterProfiler.
##
## The first argument must be the path to a file produced by one
## of the following scripts 
##   - run_edgR.R
##   - run_DESeq2.R
##
## The second argument indicates the program that was used to 
## generate the data (supported: "edgeR", "DESeq2")
##
## Authors: Lucie Khamvongsa & Jacques van Helden
## Date: August 2015
##
## Installation of the clusterProfiler package:
## source("http://bioconductor.org/biocLite.R")
## install.packages(c("httr", "qdap", "png")) ## Fix some non-declared dependencies
## biocLite()
## biocLite("clusterProfiler")
## biocLite("org.Dm.eg.db")
## biocLite("org.EcK12.eg.db")

library("clusterProfiler")

dir.main <- "~/dr-chip-rna-seq"
setwd(dir.main)

## TEMPORARY: hard-coded path. Should be adapted to read file path from the command-line arguments. 
deg.file <- ("results/rna-seq/DEG/S2_vs_WT/S2_VS_WT_sickle_se_q20_subread_featurecounts_DESeq2.tab")
deg.software <- "DESeq2"

## Load the DEG analysis result table
deg.table <- read.table(deg.file)

## Parameters
threshold <- 0.05
threshold.column <- "padj"

## Select genes according to the user-specified threshold
selected <- deg.table[,threshold.column] <= threshold
summary(selected) ## Check the number of FALSE, TRUE and NA values
selected[is.na(selected)] <- FALSE ## Replace the NA values by FALSE to avoid problems with selection

## Select up-regulated genes
summary()
