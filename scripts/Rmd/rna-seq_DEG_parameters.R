########################################################################
## Parameters required for the RNA-seq differential analysis

## In this chunk we define the project-specific parameters that are required to enable an analysis. Other (optional) parameters can be defined below to customize the thresholds, figure styles, ... 
dir.fg <- "~/fg-chip-seq/" ## Directory containing routines from France Génomique

## Load a library with the utilities
source(file.path(dir.fg, "scripts/R-scripts/deg_lib.R"))

## base directory for the report (links are defined relative to this directory)
dir.main <- "~/BeatriceRoche/"
dir.report <- "~/BeatriceRoche/reports/"
dir.base <- ".." ; ## Path to the main directory relative to the current report, which must be added to the all links
opts_knit$set(base.dir=dir.report) ## Set the working directory for knitr (generating HTML and pdf reports)
setwd(dir.report) ## Set the working directory for the console

################################################################
## Input parameters

## Sample description file (one row per sample)
sample.description.file <- "data/broche_samples.tab"

## Dsign file (pairs of samples to be compared)
design.file <- "data/broche_analyses.tab"

## Table containing the counts of reads per gene (rows) for each sample (columns) 
all.counts.table <- 'results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_all_counts.tab'

## The GTF file contains genomic annotations for the selected species
organism.name <- "Escherichia coli"
gtf.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.gtf"
gtf.source <- "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-28/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/"

gene.info.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2_gene_info.tab"

################################################################
## Output parameters
dir.DEG <-'results/DEG' ## Output directory (relative to base dir)
dir.figures <- file.path(dir.DEG, "figures")
dir.create(file.path(dir.base, dir.figures), showWarnings = FALSE, recursive = TRUE)
suffix.deg <- 'broche_project' ## Suffix for the output files
suffix.DESeq2 <- paste(sep="_", suffix.deg, "DESeq2")
suffix.edgeR <- paste(sep="_", suffix.deg, "edgeR")

## Prefix for output files concerning the whole count table (all samples together)
## prefix["general.file"] <- sub(pattern = ".tab", replacement="", all.counts.table)
prefix <- vector()
prefix["general.file"] <- file.path(dir.DEG, suffix.deg)
prefix["general.file"] <- file.path(dir.DEG, suffix.deg) ## Path prefix for the general files

################################################################
## Thresholds for differential analysis. 

thresholds <- c("padj"=0.05, ## Upper threshold on False Discovery Rate (FDR)
                "evalue"=1, ## Upper threshold on the expected number of false positives
                "FC"=1.5 ## Lower threshold on the fold-change
)

project.info <- c(
  "Customer"="Beatrice Roche",
  "Platform"="TGML/TAGC, Marseille, France",
  "Bioinfo/stat analysis"="Jacques van Helden",
  "Project start"="2015",
  "Last update"=format(Sys.time(), "%Y-%m-%d %H:%M"))
