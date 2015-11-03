########################################################################
## Parameters required for the RNA-seq differential analysis

## In this chunk we define the project-specific parameters that are required to enable an analysis. Other (optional) parameters can be defined below to customize the thresholds, figure styles, ... 
dir.fg <- "/home/lkhamvongsa/mountrsatlocal/fg-chip-seq/" ## Directory containing routines from France Génomique

## Load a library with the utilities
source(file.path(dir.fg, "scripts/R-scripts/deg_lib.R"))

## base directory for the report (links are defined relative to this directory)
dir.main <- "/home/lkhamvongsa/mountrsatlocal/Mathilde/"
dir.report <- "/home/lkhamvongsa/mountrsatlocal/Mathilde/reports/"
dir.base <- "/home/lkhamvongsa/mountrsatlocal/Mathilde"  ## Path to the main directory relative to the current report, which must be added to the all links
opts_knit$set(base.dir=dir.report) ## Set the working directory for knitr (generating HTML and pdf reports)
setwd(dir.report) ## Set the working directory for the console

################################################################
## Input parameters

## Sample description file (one row per sample)
sample.description.file <- "data/PLZF_sample_descriptions_one_raw.tab"

## Design file (pairs of samples to be compared)
design.file <- "data/PLZF_analysis_description_pairs.tab"

## Table containing the counts of reads per gene (rows) for each sample (columns) 
all.counts.table <- 'results/DEG/GMP_Mut090215_Mut130215_WT130215_WT121214.tab'

## The GTF file contains genomic annotations for the selected species
organism.name <- "Mus musculus"
gtf.file <- "genome/mm10.gtf"
gtf.source <- "ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/"

#parameter.file <- "/home/lkhamvongsa/mountrsatlocal/Mathilde/scripts/Rmd/rna-seq_DEG_parameters.R"

################################################################
## Output parameters
dir.DEG <-'results/DEG' ## Output directory (relative to base dir)
dir.figures <- file.path(dir.DEG, "figures")
dir.create(file.path(dir.base, dir.figures), showWarnings = FALSE, recursive = TRUE)
suffix.deg <- 'mathilde_project' ## Suffix for the output files
suffix.DESeq2 <- paste(sep="_", suffix.deg, "DESeq2")
suffix.edgeR <- paste(sep="_", suffix.deg, "edgeR")

## Prefix for output files concerning the whole count table (all samples together)
## prefix["general.file"] <- sub(pattern = ".tab", replacement="", all.counts.table)
prefix <- vector()
prefix["general.file"] <- file.path(dir.DEG, suffix.deg)
prefix["general.file"] <- file.path(dir.DEG, suffix.deg) ## Path prefix for the general files

################################################################
## Thresholds for differential analysis. 
#"evalue"=2, ## Upper threshold on the expected number of false positives
thresholds <- c("padj"=0.05, ## Upper threshold on False Discovery Rate (FDR)
                "FC"=1.5 ## Lower threshold on the fold-change
)


project.info <- c(
  "Customer"="Mathilde CRCM",
  "Platform"="TGML/TAGC, Marseille, France",
  "Bioinfo/stat analysis"="Mathilde Poplineau& Lucie Khamvongsa",
  "Project start"="2015",
  "Last update"=Sys.Date())
