################################################################
## R code to detect differentially expressed genes with a combination 
## of several  BioConductor libraries, including edgeR and DESeq2.
##
## Authors: Jeanne Chèneby & Justine Long (Initial script), revised and completed by Jacques van Helden
## First version: 2015-07
## Last update: 2015-08-23
################################################################


##Load libraries
library("edgeR", warn.conflicts = FALSE, quietly=TRUE)
library("DESeq2", warn.conflicts = FALSE, quietly=TRUE, verbose=FALSE)
library("limma", warn.conflicts = FALSE, quietly=TRUE) ## Required for vennCounts and vennDiagram
library(gplots, warn.conflicts = FALSE, quietly=TRUE) ## Required for heatmaps.2
library(RColorBrewer, warn.conflicts = FALSE, quietly=TRUE)
library("GenomicFeatures")
library("clusterProfiler")
library("stats4bioinfo")



## Load library of functions for differential analysis
dir.fg <- "~/fg-chip-seq/"
#dir.fg <- "~/mountrsatlocal/fg-chip-seq/"
source(file.path(dir.fg, "scripts/R-scripts/deg_lib.R"))

## Define parameters
run.param <- list()
run.param$exploratory.plots <- FALSE
deg.tools <- c("edgeR", "DESeq2")

verbosity <- 1
save.image <- TRUE ## Save memory image in an RData file

## The only argument is the file containing all the parameters for the analysis
if (is.na(commandArgs(trailingOnly = FALSE)[6])) {
  #   stop("Parameter file must be passed as command-line argument.")
  ## Elements that should be added to the parameters
  #org <- "dm"
  #org <- "dm_sepsis"
  #org <- "eco"
  if (org == "eco") {
    ## TEMPORARY FOR DEBUGGING: 
    source ("~/BeatriceRoche/config_files/roche-loiseau_config.R")
#     dir.main <- "~/BeatriceRoche/"
#     setwd(dir.main)
#     r.params.path <- "results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_params.R"  
#     org.db <- "org.EcK12.eg.db" ## Should be added to parameters
#     gene.info.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2_gene_info.tab"
#     organism.names <- c("name" = "Escherichia coli",
#                         "clusterProfiler" = NA,
#                         "kegg"="eco")
#     gtf.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.gtf"
#     gtf.source <- "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-28/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/"
#     #   pet.gene <- "b2531"
#     genes.of.interest <- c("b2531")
#     go.map.file <- 'genome/Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2_gene_GO.tab'
#     go.description.file <- "genome/GO_description.tab"
  } else if (org=="dm") {
    dir.main <- "~/dr-chip-rna-seq/"
    setwd(dir.main)
    r.params.path <- "results/rna-seq/DEG/sickle_se_q20_subread_featurecounts_params.R"  
    sample.description.file <- "data/rna-seq/abdA_overexpr_sample_descriptions.tab" 
    design.file <- "data/rna-seq/analysis_description.tab"
    gtf.file <- "genomes/ftp.ensemblgenomes.org/pub/metazoa/release-28/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.gtf"
    gtf.source <- "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-28/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.gtf"
    dir.DEG <- "results/rna-seq/DEG/"
    suffix.deg <- "sickle_se_q20_subread_featurecounts"
    suffix.DESeq2 <- paste(sep="_", suffix.deg, "DESeq2")
    suffix.edgeR <- paste(sep="_", suffix.deg, "edgeR")
    org.db <- "org.Dm.eg.db"
    organism.names <- c("name" = "Drosophila melanogaster",
                        "clusterProfiler" = "fly",
                        "kegg"="dme")
  }else if (org=="dm_sepsis") {
    dir.main <- "/home/lkhamvongsa/mountrsatlocal/droso_sepsis"
    setwd(dir.main)
    r.params.path <- "results/DEG/sickle-se-q20_subread_featurecounts_params.R"  
    sample.description.file <- "data/Trauma_induced_sepsis_sample_description.tab"
    design.file <- "data/Trauma_induced_sepsis_analysis_description.tab"
    gtf.file <- "genome/genes.gtf"
    gtf.source <- "/home/lkhamvongsa/mountrsatlocal/workspace/genomes/Drosophila_melanogaster_UCSC_dm3/dm3/Annotation/Genes/genes.gtf"
    dir.DEG <- "results/DEG/"
    suffix.deg <- "sickle-se-q20_subread_featurecounts"
    suffix.DESeq2 <- paste(sep="_", suffix.deg, "DESeq2")
    suffix.edgeR <- paste(sep="_", suffix.deg, "edgeR")
    org.db <- "org.Dm.eg.db"
    organism.names <- c("name" = "Drosophila melanogaster",
                        "clusterProfiler" = "fly",
                        "kegg"="dme")
  }else {
    message(paste("No info for org", org))
  }
} else {
  r.params.path <- commandArgs(trailingOnly = FALSE)[6]  
}
verbose(paste("Reading parameters", r.params.path), 1)
source(r.params.path)
#setwd(dir.main)
# getwd()

################################################################
## Check parameters + define additional ones

## Experimental: we would like to compare gene lists obtained by applying thresholds on different scores. 
## At the end we can choose any combination of these threshold to select the relevant gene list.
## Choose default value for thresholds if they were not defined in the config file
if (!exists("thresholds")) {
  thresholds <- c("padj"=0.05, "FC"=1.5) 
}


## Define a color palette for heatmaps. I like this Red-Blue palette because 
## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
## - it can be seen by people suffering from red–green color blindness.
cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))

## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
epsilon <- 0.01

## Read the sample description file, which indicates the 
## condition associated to each sample ID.
verbose("Reading sample descriptions", 1)

sample.desc <- read.delim(sample.description.file, sep="\t", 
                          comment=";", header = TRUE, row.names=1)
sample.ids <- row.names(sample.desc)

## Experimental conditions
sample.conditions <- as.vector(sample.desc[,1]) ## Condition associated to each sample
names(sample.conditions) <- sample.ids
# print(sample.conditions)

## Build sample labels by concatenating their ID and condition
sample.desc$label <- paste(sep="_", sample.ids, sample.conditions)

## Define a specific color for each distinct condition
conditions <- unique(sample.conditions) ## Set of distinct conditions
cols.conditions <- brewer.pal(max(3, length(conditions)),"Dark2")[1:length(conditions)]
names(cols.conditions) <- conditions
# print(cols.conditions)

## Define a color per sample according to its condition
sample.desc$color <- cols.conditions[sample.conditions]
# names(cols.samples) <- sample.ids
# print(cols.samples)

## Read the design file, which indicates the anlayses to be done.
## Each row specifies one differential expression analysis, which 
## consists in comparing two conditions. 
verbose("Reading design", 1)
design <- read.delim(design.file, sep="\t", 
                     comment=";", header = T, row.names=NULL)
# print(design)

## Prefix for output files concerning the whole count table (all samples together)
## prefix["general.file"] <- sub(pattern = ".tab", replacement="", all.counts.table)
prefix <- vector()
prefix["general.file"] <- file.path(dir.DEG, suffix.deg)
# print(prefix)

## Read the count table
verbose("Loading count table", 1)
all.counts <- read.delim(all.counts.table, row.names=1, sep="\t")
# names(all.counts)
# dim(all.counts)


## Check that the header of all.counts match the sample IDs
ids.not.found <- setdiff(sample.ids, names(all.counts)) ## Identify sample IDs with no column in the count table
if (length(ids.not.found) > 0) {
  stop(paste(length(ids.not.found), "Missing columns in count table", all.counts.table), paste(sep="; ", ids.not.found))
}

verbose("Computing count-derived metrics (log-transformed)", 1)

## The next step is useful only with htseq-count results, 
## in order to suppress the statistics that are inserted within the count files themselves.
## Statistics on reads that could not be mapped to genes for different reasons (intergenic, ambiguous, not unique, ...)
all.counts.htseq.stats <- all.counts[grep(pattern = "^__", x = row.names(all.counts)), ]
# dim(all.counts.htseq.stats)
all.counts.mapped <- all.counts[grep(invert=TRUE, pattern = "^__", x = row.names(all.counts)), ]
# dim(all.counts.mapped)

## Add an epsilon to 0 values only, in order to enable log-transform and display on logarithmic axes
all.counts.mapped.epsilon <- all.counts.mapped
all.counts.mapped.epsilon[all.counts.mapped==0] <- epsilon

## Log-transformed data for some plots. 
all.counts.mapped.log10 <- log10(all.counts.mapped.epsilon)


################################################################
## Load gene information from the GTF file 
## (should be the same as used to count tags per gene)
if (exists("gtf.file")) {
  
  ## Jacques or Lucie, TO DO : check if the method makeTxDbFromGFF allows to get gene names and descriptions.
  ## In principle it should be possible since this info is in the GTF file. 
  ## If not, we migh rewite a parser for GTF. 
  
  verbose(paste("Loading gene information from GTF file", gtf.file))
  #   library(rtracklayer)  # for the import() function
  #   gr <- import(gtf.file)
  txdb <- makeTxDbFromGFF(file=gtf.file,
                          organism=organism.names["name"],
                          dataSource=gtf.source)
  #seqlevels(txdb)
  #   transcripts <- transcriptsBy(txdb, by="gene")
  #   transcript.table <- as.data.frame(transcripts@unlistData)
  #   cds <- cdsBy(txdb, by="gene")
  #   cds.table <- as.data.frame(cds@unlistData)
  #   exons <- exonsBy(txdb, by="gene")
  #   exon.table <- as.data.frame(exons@unlistData)
  all.genes <- genes(txdb)
  gene.info <- as.data.frame(all.genes)
  gene.info$name <- gene.info$gene_id
  gene.info$entrez.id <- NA
  gene.info$description <- "no description"
} else {
  verbose(paste("No GTF file has been specified"))
  g <- nrow(all.counts.mapped)
  gene.info <- data.frame("seqnames"=rep(NA, times=g),
                          "start"=rep(NA, times=g),
                          "end"=rep(NA, times=g),
                          "width"=rep(NA, times=g),
                          "strand"=rep(NA, times=g),
                          "gene_id"=row.names(all.counts.mapped),
                          "name"=row.names(all.counts.mapped),
                          "entrez.id" = rep(NA, times=g),
                          "description"=rep("no description", times=g))
  row.names(gene.info) <- row.names(all.counts.mapped)
}
# View(gene.info)



################################################################
## TEMPORARY: load additional gene information from a tab-delimited file generated with RSAT, 
## because for bacterial genomes there is no obvious way to obtain EntrezIDs.
if (exists("gene.info.file")) {
  verbose(paste("Loading gene information from tab-delimited file", gene.info.file))
  
  gene.info.rsat <- read.delim(gene.info.file, sep="\t", skip=3, header=1)  
  ## A bit tricky: we need to discard rows starting with ";", but we cannot use ";" as comment.char since it is found in descriptions and as separator for names
  comment.lines <- grep(gene.info.rsat[,1], pattern = ";")
  gene.info.rsat <- gene.info.rsat[-comment.lines,]
  row.names(gene.info.rsat) <- gene.info.rsat[,1]
  gene.info.rsat <- gene.info.rsat[row.names(gene.info), ] ## Ensure gene.info.rsat contains same IDs in the same order as gene.info
  
  
  ## Extract additional attributes for the gene.info table
  gene.info$name <- gene.info.rsat[row.names(gene.info), "name"]
  gene.info$entrez.id <- gene.info.rsat[row.names(gene.info), 2]
  gene.info$description <- gene.info.rsat[row.names(gene.info), "description"]
  gene.info$names <- gene.info.rsat[row.names(gene.info), "names"]
  # dim(gene.info.rsat)
  # names(gene.info.rsat)
  # View(gene.info.rsat)
} 
# View(gene.info)

## Export the gene information table to keep a trace of what has been used. 
gene.info.out <- paste(sep="", prefix["general.file"], "_gene_descriptions.tab")
verbose(paste("Exporting gene information table", gene.info.out), 1)
write.table(x = gene.info, row.names = FALSE,
            file = gene.info.out, sep = "\t", quote=FALSE)
verbose(paste(sep="", "\tGene info table\t", gene.info.out), 1)


################################################################
## Compute sample-wise statistics on mapped counts
################################################################
verbose("Computing statistics per sample", 1)

stats.per.sample <- cbind(  
  sample.desc[names(all.counts.mapped), ],
  data.frame(
    "sum" = apply(all.counts.mapped, 2, sum, na.rm=TRUE),
    "mean" = apply(all.counts.mapped, 2, mean, na.rm=TRUE),
    "min" = apply(all.counts.mapped, 2, min, na.rm=TRUE),
    "perc05" = apply(all.counts.mapped, 2, quantile, probs=0.05, na.rm=TRUE),
    "perc25" = apply(all.counts.mapped, 2, quantile, probs=0.25, na.rm=TRUE),
    "median" = apply(all.counts.mapped, 2, median, na.rm=TRUE),
    "perc75" = apply(all.counts.mapped, 2, quantile, probs=0.75, na.rm=TRUE),
    "perc95" = apply(all.counts.mapped, 2, quantile, probs=0.95, na.rm=TRUE),
    "max" = apply(all.counts.mapped, 2, max, na.rm=TRUE)
  )
)
stats.per.sample$Mreads <- round(stats.per.sample$sum/1e6, digits = 1)
  
## Count number and the fraction of samples with counts below the mean. 
## This shows the impact of very large counts: in my test samples, 
## 85% of the samples have a value below the mean (i.e. the mean is at the percentile 85 !)
stats.per.sample$below.mean <- apply(t(all.counts.mapped) < stats.per.sample$mean, 1, sum, na.rm=TRUE)
stats.per.sample$fract.below.mean <- stats.per.sample$below.mean/nrow(all.counts.mapped)
# View(stats.per.sample)




################################################################
## Compute the counts per million reads 
################################################################
verbose("Computing CPMs", 1)

## Note: the default normalization criterion (scaling by libbrary sum) 
## is questionable because it is stronly sensitive to outliers 
## (very highly expressed genes).  A more robust normalisation criterion 
## is to use the 75th percentile, or the median. We use the median, somewhat arbitrarily, 
## beause it gives a nice alignment on the boxplots.
cpms.libsum <- cpm(all.counts.mapped.epsilon)    ## Counts per million reads, normalised by library sum
cpms.perc75 <- cpm(all.counts.mapped.epsilon, lib.size = stats.per.sample$perc75)    ## Counts per million reads, normalised by 75th percentile
cpms.perc95 <- cpm(all.counts.mapped.epsilon, lib.size = stats.per.sample$perc95)    ## Counts per million reads, normalised by 95th percentile
cpms.median <- cpm(all.counts.mapped.epsilon, lib.size = stats.per.sample$median)    ## Counts per million reads, normalised by sample-wise median count
#cpms <- cpms.median ## Choose one normalization factor for the CPMs used below
cpms <- cpms.perc75 ## Choose one normalization factor for the CPMs used below
cpms.log10 <- log10(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts
cpms.log2 <- log2(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts

stats.per.sample$cpm.mean <- apply(cpms, 2, mean)
stats.per.sample$log2.cpm.mean <- apply(cpms.log2, 2, mean)
stats.per.sample$log10.cpm.mean <- apply(cpms.log10, 2, mean)

################################################################
## Export stats per sample
#
# names(stats.per.sample)
# head(stats.per.sample)
verbose("Exporting stats per sample", 1)
sample.summary.file <- paste(sep="", prefix["general.file"], "_summary_per_sample.tab")
write.table(x = stats.per.sample, row.names = TRUE, col.names=NA, 
            file = sample.summary.file, sep = "\t", quote=FALSE)
verbose(paste(sep="", "\tSummary per sample\t", sample.summary.file), 1)

################################################################
## Draw some generic plots to display sample-wise statistics.
################################################################

plot.files <- sample.description.plots(
  sample.desc, stats.per.sample, dir.DEG, 
  exploratory.plots=run.param$exploratory.plots)

################################################################
## Analyse between-replicate reproducibility
################################################################
verbose("Plotting betwen-replicate comparisons", 1)
#conditions.with.replicate <- setdiff(conditions, "CI_CI_0H")
cond <- conditions[1] ## Choose first condition for testing without the loop
for (cond in conditions) {
#for (cond in conditions.with.replicate) {
  verbose(paste(sep="", "\tcondition\t", cond), 2)
  
  max.rep.to.plot <- 5 ## Restrict the number of replicates to plot, for the sale of readability
  
  ## Create a specific result directory for this condition
  dir.condition <- file.path(dir.DEG, "per_condition", cond)
  dir.create(path = dir.condition, showWarnings = FALSE, recursive = TRUE)
  
  ## Select the specific data for the current condition (samples, counts, CPMs)
  #all.counts.mapped.with.replicate <- setdiff(names(all.counts.mapped), "CI_CI_0H_1")
  current.samples <- names(all.counts.mapped)[sample.conditions == cond]
  nrep <- length(current.samples)
  
  current.counts <- all.counts.mapped[,current.samples]
  current.counts[current.counts==0] <- epsilon
  current.counts.mean <- apply(current.counts, 1, mean)
  current.counts.var <- apply(current.counts, 1, var)
  current.counts.var[current.counts.var==0] <- min(current.counts.var[current.counts.var>0])/100
  
  current.cpms <- cpms[,current.samples] + epsilon
  current.cpms[current.cpms==0] <- epsilon
  current.cpm.mean <- apply(current.cpms, 1, mean)
  current.cpm.var <- apply(current.cpms, 1, var)
  current.cpm.var[current.cpm.var==0] <- min(current.cpm.var[current.cpm.var>0])/100
  
  ################################################################
  ## Plot pairwise comparisons between replicates (result file can be heavy to open)
  pdf(file= file.path(dir.condition, paste(sep = "", "between-replicate_counts_plot_", cond, ".pdf")), width=10, height=10)
  # dim(current.counts[,1:min(max.rep.to.plot, nrep)])
  plot(current.counts[,1:min(max.rep.to.plot, nrep)], 
       log="xy", col=cols.conditions[cond], 
       #       panel.first=grid(col="#BBBBBB", lty="solid"), ## this does not work, I should check why
       main=paste(cond, " ; raw counts per replicate (log scale)")
  )
  silence <- dev.off()
  
  ################################################################
  ## Plot mean versus variance of CPMs for the current condition.
  ## BEWARE! This is the between-replicate variance computed for each gene separately, 
  ## which is quite different from the smoothed estimate of variance used by 
  ## DESeq2 or edgeR for the negative binomial.
  ## 
  ## This "real" variance plot  can be useful to highlight some genes that show 
  ## an extremely high variance compared to other genes. 
  ##
  ## We draw the plot in linear + log scales. Linear scales better highlight the outliers, 
  ## log scale better shows the general distribution, which covers several orders of magnitude.
  log.axes <- "xy"
  for (log.axes in c("", "xy")) {
    pdf(file= file.path(dir.condition, paste(sep = "", "CPM_variance-mean_plot_", cond, "_",  log.axes, ".pdf")), width=10, height=10)
    plot(current.cpm.mean, 
         current.cpm.var, 
         log=log.axes, col=cols.conditions[cond], 
         panel.first=grid(lty="solid", col="#DDDDDD"),
         xlab=paste("Mean CPM per gene for condition", cond),
         ylab=paste("Between replicate CPM variance per gene for condition", cond),
         main=paste(cond, " ; CPM variance/Mean plot"))
    abline(a=0,b=1, lty="dashed", col="green", lwd=2) ## Milestone for Poisson distributions: var = mean
    silence <- dev.off() 
  }
}

################################################################
## Run differential expression analysis
################################################################
verbose("Starting differential analysis", 1)

## Iterate over analyses
i <- 1
comparison.results <- design
comparison.results$prefixes <- paste(sep="_", design$cond1, "vs", design$cond2)
for (i in 1:nrow(design)) {
  
  ## Identify samples for the first condition
  cond1 <- as.vector(design[i,1])  ## First condition for the current comparison
  samples1 <- sample.ids[sample.conditions == cond1]
   if (length(samples1) < 2) {
     stop(paste("Cannot perform differential analysis. The count table contains less than 2 samples for condition", cond1))
   }
  
  ## Identify samples for the second condition
  cond2 <- as.vector(design[i,2])  ## Second condition for the current comparison
  samples2 <- sample.ids[sample.conditions == cond2]
   if (length(samples2) < 2) {
     stop(paste("Cannot perform differential analysis. The count table contains less than 2 samples for condition", cond2))
   }
  
  verbose(paste(sep="", "\tDifferential analysis\t", i , "/", nrow(design), "\t", cond1, " vs ", cond2), 1)
  
  ## Create a specific result directory for this differential analysis
  prefix["comparison"] <- comparison.results$prefixes[i]
  dir.analysis <- file.path(dir.DEG, paste(sep="", prefix["comparison"]))
  comparison.results$dir.analysis <- dir.analysis
  dir.create(path = dir.analysis, showWarnings = FALSE, recursive = TRUE)
  dir.figures <- file.path(dir.analysis, "figures")
  comparison.results$dir.figures <- dir.figures
  dir.create(path = dir.figures, showWarnings = FALSE, recursive = TRUE)
  prefix["comparison_file"] <- file.path(dir.analysis, prefix["comparison"])
  prefix["comparison_figure"] <- file.path(
    dir.figures, 
    paste(sep="", prefix["comparison"], "_",  suffix.deg))
  
  ## Select counts for the samples belonging to the two conditions
  current.samples <- c(samples1, samples2)
  current.counts <- all.counts.mapped[,current.samples]
  # dim(current.counts)  ## For test
  # names(current.counts)
  
  if (sum(!names(current.counts) %in% sample.ids) > 0) {
    stop("Count table contains column names without ID in sample description file.")
  }
  
  ## Define conditions and labels for the samples of the current analysis
  current.sample.conditions <- sample.conditions[current.samples]
  current.labels <- paste(current.sample.conditions, names(current.counts), sep="_")
  
  ## Initiate a result table with the CPMs and derived statistics
  all.gene.ids <- row.names(cpms)
  result.table <- data.frame("gene_id" = all.gene.ids,
                             "name"=gene.info[all.gene.ids,"name"])
  row.names(result.table) <- all.gene.ids
  result.table$entrez.id <- gene.info[all.gene.ids,"entrez.id"]
  result.table$description <- gene.info[all.gene.ids,"description"]
  
  result.table <- cbind(result.table, counts=current.counts) ## Include the original counts in the big result table
  result.table <- cbind(result.table, cpm=current.cpms) ## Include CPMs in the big result table
  
  ## Tag genes detected in less than min.rep samples, which is defined as 
  ## the minimal number of replicates per condition.
  min.rep <- min(length(samples1), length(samples2))
  result.table$undetected <- rowSums(current.counts > 1) < min.rep
  # table(result.table$undetected)
  # dim(current.counts)
  
  result.table$cpm.mean <- apply(cpms,1, mean)
  result.table$cpm1.mean <- apply(as.data.frame(cpms[,samples1]),1, mean)
  result.table$cpm2.mean <- apply(as.data.frame(cpms[,samples2]),1, mean)
  result.table$A = log2(result.table$cpm1.mean*result.table$cpm2.mean)/2
  result.table$M = log2(result.table$cpm1.mean/result.table$cpm2.mean)
  result.table$cpm.median <- apply(cpms,1, median)
  result.table$cpm1.median <- apply(as.data.frame(cpms[,samples1]),1, median)
  result.table$cpm2.median <- apply(as.data.frame(cpms[,samples2]),1, median)
  result.table$cpm.min <-  apply(cpms,1, min)
  result.table$cpm1.min <- apply(as.data.frame(cpms[,samples1]),1, min)
  result.table$cpm2.min <- apply(as.data.frame(cpms[,samples2]),1, min)
  result.table$cpm.max <-  apply(cpms,1, max)
  result.table$cpm1.max <- apply(as.data.frame(cpms[,samples1]),1, max)
  result.table$cpm2.max <- apply(as.data.frame(cpms[,samples2]),1, max)
  result.table$cpm.sd <-  apply(cpms,1, sd)
  result.table$cpm1.sd <- apply(as.data.frame(cpms[,samples1]),1, sd)
  result.table$cpm2.sd <- apply(as.data.frame(cpms[,samples2]),1, sd)
  result.table$cpm.var <-  apply(cpms,1, var)
  result.table$cpm1.var <- apply(as.data.frame(cpms[,samples1]),1, sd)
  result.table$cpm2.var <- apply(as.data.frame(cpms[,samples2]),1, sd)
  # View(result.table)
  
  ################################################################
  ## DESeq2 analysis
  ################################################################
  
  verbose("\t\tDESeq2 analysis", 2)
  
  ## Path prefix to save DESeq2 result files
  prefix["DESeq2_file"] <- paste(sep="", prefix["comparison_file"], "_", suffix.DESeq2)
  prefix["DESeq2_figure"] <- paste(sep="", prefix["comparison_figure"], "_", suffix.DESeq2)
  
  ## Create a DESeqDataSet object from the count table + conditions
  condition <- as.factor(as.vector(current.sample.conditions))
  deseq2.dds <- DESeqDataSetFromMatrix(
    countData = current.counts, 
    colData = DataFrame(condition),
    ~ condition)
  
  
  ## Indicate that second condition is the reference condition. 
  ## If not done, the conditions are considered by alphabetical order, 
  ## which may be misleading to interpret the log2 fold changes. 
  deseq2.dds$condition <- relevel(deseq2.dds$condition, ref=cond2) 
  
  ## Run the differential analysis
  deseq2.dds <- DESeq(deseq2.dds)      ## Differential analysis with negbin distrib
  deseq2.res <- results(deseq2.dds, independentFiltering=FALSE, pAdjustMethod = "BH")  ## Collect the result table
  
  deseq2.result.table <- data.frame(
    "gene.id" = row.names(deseq2.res),
    "mean" = deseq2.res$baseMean,
    "log2FC" = deseq2.res$log2FoldChange,
    "pvalue" = deseq2.res$pvalue,
    "padj" = deseq2.res$padj)
  deseq2.result.table <- complete.deg.table(
    deseq2.result.table, 
    paste(sep="_", "DESeq2", prefix["comparison"]),
    sort.column = "padj",
    thresholds=thresholds,
    round.digits = 3,
    dir.figures=dir.figures)
  result.table <- cbind(result.table, "DESeq2" = deseq2.result.table[row.names(result.table),])
  # dim(deseq2.result.table)
  # dim(deseq2.result.table)
  # names(deseq2.result.table)
  # View(deseq2.result.table)
  # View(result.table)
  
  ## Save the completed DESeq2 result table
  deseq2.result.file <- paste(sep = "", prefix["DESeq2_file"], ".tab")
  comparison.results[i,"deseq2"] <- deseq2.result.file
  write.table(x = deseq2.result.table, 
              row.names = FALSE, file = deseq2.result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tDESeq2 result file\t", deseq2.result.file), 1)
  
  ################################################################
  ## Export DESeq2 plots
  verbose("\t\tExporting DESeq2 plots", 2)
  
  
  # Transform the raw discretely distributed counts to apply 
  # some multivariate methods such as PCA, clustering etc.
  deseq2.rld <- rlogTransformation(deseq2.dds, blind=TRUE)
  
  ## MA plot
  verbose("\t\t\tDESeq2 MA plot", 2)
  MA.ymax <- min(4,max(na.omit(deseq2.res$log2FoldChange)))
  MA.ymin <- max(-4,min(na.omit(deseq2.res$log2FoldChange)))
  MA.ylim <- c(MA.ymin, MA.ymax)
  pdf(file=paste(sep = "", prefix["DESeq2_figure"], "_MAplot.pdf"))
  plotMA(deseq2.dds, main=paste("DESeq2", prefix["comparison"]), ylim=MA.ylim)
  grid(lty="solid", col="#CCCCCC")
  quiet <- dev.off()
  
  # PCA plot of the samples
  verbose("\t\t\tDESeq2 PCA plot", 2)
  pdf(file=paste(sep = "", prefix["DESeq2_figure"], "_PCAplot.pdf"))
  print(plotPCA(deseq2.rld, intgroup=c("condition")))
  quiet <- dev.off()
  
  ################################################################
  ## edgeR analysis
  ################################################################
  verbose("\t\tedgeR analysis", 2)
  
  prefix["edgeR_file"] <- paste(sep="", prefix["comparison_file"], "_", suffix.edgeR)
  prefix["edgeR_figure"] <- paste(sep="", prefix["comparison_figure"], "_", suffix.edgeR)
  
  ## Convert the count table in a DGEList structure and compute its parameters.
  d <- DGEList(counts=current.counts, group=sample.conditions[names(current.counts)])
  d$samples$group <- relevel(d$samples$group, ref=cond2) ## Ensure that condition 2 is considered as the reference
  d <- calcNormFactors(d, method="RLE")                 ## Compute normalizing factors
  d <- estimateCommonDisp(d, verbose=FALSE)             ## Estimate common dispersion
  d <- estimateTagwiseDisp(d, verbose=FALSE)            ## Estimate tagwise dispersion
  
  ################################################################
  ## Detect differentially expressed genes by applying the exact 
  ## negative binomial test from edgeR package.
  edger.de <- exactTest(d, pair=c(cond2, cond1))      ## Run the exact negative binomial test
  
  ## Sort genes by increasing p-values, i.e. by decreasing significance
  edger.tt <- topTags(edger.de, n=nrow(d), sort.by = "PValue")
  
  ## Complete the analysis of edgeR result table
  edger.result.table <- data.frame("gene.id" = row.names(edger.tt$table),
                                   "mean"=edger.tt$table$logCPM,
                                   "log2FC"=edger.tt$table$logFC,
                                   "pvalue"=edger.tt$table$PValue,
                                   "padj"=edger.tt$table$FDR)
  edger.result.table <- complete.deg.table(
    deg.table = edger.result.table, 
    table.name = paste(sep="_","edgeR",prefix["comparison"]),
    sort.column = "padj",
    thresholds=c("padj"=0.05, "evalue"=1, "FC"=1.5),
    round.digits = 3,
    dir.figures=dir.figures)
  result.table <- cbind(result.table, "edgeR" = edger.result.table[row.names(result.table),])
  # dim(edger.result.table)
  # dim(edger.result.table)
  # names(edger.result.table)
  # View(edger.result.table)
  
  ## Export edgeR result table
  edger.result.file <- paste(sep="", prefix["edgeR_file"], ".tab")
  comparison.results[i,"edger"] <- edger.result.file
  
  write.table(x = edger.result.table, 
              row.names = FALSE, file = edger.result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tedgeR result file\t", edger.result.file), 1)
  
  ################################################################
  ## Export edgeR plots
  
  verbose("\t\tGenerating edgeR plots", 2)
  
  ## Smear plot
  verbose("\t\t\tedgeR Smear plot", 3)
  edger.deg <- rownames(edger.tt$table)[edger.tt$table$FDR < thresholds["padj"]] ## List of differentially expressed genes reported by edgeR
  pdf(file=paste(sep="", prefix["edgeR_figure"], "_plotSmear.pdf"))
  plotSmear(d, de.tags = edger.deg)
  quiet <- dev.off()
  
  ## MDS plot (Multidimensional scaling plot of distances between gene expression profiles)
  verbose("\t\t\tedgeR MDS plot", 3)
  pdf(file=paste(sep="", prefix["edgeR_figure"], "_plotMDS.pdf"))
  #  pdf(file=file.path(dir.figures, paste(sep="", "edgeR_MDS_plot_", prefix["comparison"], ".pdf")))
  plotMDS(d, labels=current.labels, 
          col=c("darkgreen","blue")[factor(sample.conditions[names(current.counts)])])
  quiet <- dev.off()
  
  # Mean-variance relationship
  verbose("\t\t\tedgeR mean-variance relationship", 3)
  pdf(file=paste(sep="", prefix["edgeR_figure"], "_plotMeanVar.pdf"))
  plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
  quiet <- dev.off()
  
  
  # BCV (Biological Coefficient of Variation)
  verbose("\t\t\tedgeR BCV plot", 3)
  pdf(file=paste(sep="", prefix["edgeR_figure"], "_plotBCV.pdf"))
  #pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotBCV_", prefix["comparison"], ".pdf")))
  plotBCV(d)
  quiet <- dev.off()
  
  ################################################################
  ## Compare gene selections by edgeR and DESeq2
  verbose("\t\tComparing DESeq2 and edgeR results", 2)    
  ## Add tags indicating the combinations of edgeR/DESeq2 selections for the DEG, and then for each threshold criterion taken separately.
  selection.columns <- vector()
  for (s in c("DEG", names(thresholds))) {
    if (s == "DEG") {
      selection.column <- "DEG"
    } else {
      selection.column <- paste(sep="_", s, thresholds[s])
    }
    selection.columns[s] <- selection.column
    
    ## Tag gnes selected by both edgeR and DESeq2
    result.table[, paste(sep="_", s, "edgeR_and_DESeq2")] <- 
      1*(result.table[,paste(sep="", "edgeR.", selection.column)] 
         & result.table[,paste(sep="", "DESeq2.", selection.column)])
    
    ## Tag genes selected by edgeR but not DESeq2
    result.table[, paste(sep="_", s, "edgeR_not_DESeq2")] <- 
      1*(result.table[,paste(sep="", "edgeR.", selection.column)] 
         & !result.table[,paste(sep="", "DESeq2.", selection.column)])
    
    ## Tag genes selected by DESeq2 but not edgeR
    result.table[, paste(sep="_", s, "DESeq2_not_edgeR")] <- 
      1*(!result.table[,paste(sep="", "edgeR.", selection.column)] 
         & result.table[,paste(sep="", "DESeq2.", selection.column)])
    
    ## Tag genes selected by DESeq2 but not edgeR
    result.table[, paste(sep="_", s, "none")] <- 
      1*(!result.table[,paste(sep="", "edgeR.", selection.column)] 
         & !result.table[,paste(sep="", "DESeq2.", selection.column)])
  }
  # View(result.table)
  
  ## Save result table
  result.file <- paste(sep = "", 
                       prefix["comparison_file"], 
                       "_", suffix.deg, "_DESeq2_and_edgeR.tab")
  comparison.results[i,"result.table"] <- result.file
  
  write.table(x = result.table, row.names = FALSE,
              file = result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tMerged result file\t", result.file), 1)
  
  
  ################################################################
  ## Plots to compare edgeR and DESeq2 results
  
  ## Define point colors according to the test results
  both <- result.table$DEG_edgeR_and_DESeq2 == 1
  edgeR.only <- result.table$DEG_edgeR_not_DESeq2 == 1
  DESeq2.only <- result.table$DEG_DESeq2_not_edgeR == 1
  none <- result.table$DEG_none
  gene.palette <- c("both"="darkgreen", 
                    "none"="#BBBBBB",
                    "edgeR.only" = "red",
                    "DESeq2.only" = "orange")
  gene.colors <- rep(x = gene.palette["none"], times = nrow(result.table))
  gene.colors[both] <- gene.palette["both"]
  gene.colors[edgeR.only] <-  gene.palette["edgeR.only"]
  gene.colors[DESeq2.only] <-  gene.palette["DESeq2.only"]
  
  ################################################################
  ## P-value comparison plot
  ## Compare DESeq2 and edgeR nominal p-values
  #png(file=paste(sep = "", prefix["comparison_figure"], suffix.deg, "_DESeq2_vs_edgeR_pvalues", ".png"))
  pdf(file=paste(sep = "", prefix["comparison_figure"], suffix.deg, "_DESeq2_vs_edgeR_pvalues", ".pdf"))
  plot.cex=0.7
  plot(result.table$edgeR.pvalue, 
       result.table$DESeq2.pvalue,
       pch=20, 
       cex=plot.cex,
       xlab="edgeR nominal p-value (log scale)", ylab="DESeq2 nominal p-value (log scale)",
       log="xy", main=paste(sep=" ", cond1, "vs", cond2, "; P-value comparisons"),
       col=gene.colors,
       panel.first=grid(lty="solid", col="#DDDDDD"))
  abline(a=0, b=1)
  legend("topleft", 
         pch=1, cex=1, 
         bg="white", bty="o",
         col=gene.palette[c("both", "DESeq2.only", "edgeR.only", "none")], 
         legend = c(paste(sum(both, na.rm=TRUE), "both"),
                    paste(sum(DESeq2.only, na.rm=TRUE), "DESeq2 only"),
                    paste(sum(edgeR.only, na.rm=TRUE), "edgeR only"),
                    paste(sum(none, na.rm=TRUE), "none")))
  silence <- dev.off()
  
  ################################################################
  ## Draw Venn diagram with number of genes declared significant with edgeR and DESeq2, resp
  
  ## Lists of differentially expressed genes (DEG) based on all thresholds together
  venn.counts.deg <- vennCounts(result.table[,c("edgeR.DEG", "DESeq2.DEG")])
  
  pdf(file=paste(sep = "", prefix["comparison_figure"], suffix.deg, "_DESeq2_vs_edgeR_Venn_DEG", ".pdf"))
  #   pdf(file= file.path(
  #     dir.figures, 
  #     paste(sep = "", "DESeq2_vs_edgeR_Venn",
  #           "_DEG", 
  #           "_", prefix["comparison"], ".pdf")))
  vennDiagram(venn.counts.deg, cex=1, 
              main=paste(sep='; ', 
                         paste(sep=" ", cond1, "vs", cond2, "DEG")))
  silence <- dev.off()
  
  
  ## For each threshold separately
  for (s in names(thresholds)) {   
    venn.counts.one.threshold <- vennCounts(
      result.table[,c(paste(sep="", "edgeR.",s,"_", thresholds[s]), 
                      paste(sep="", "DESeq2.",s,"_", thresholds[s]))])
    
    pdf(file=paste(sep = "", prefix["comparison_figure"], suffix.deg, "_DESeq2_vs_edgeR_Venn", 
                   "_", s, "_", thresholds[s], ".pdf"))
    #     pdf(file= file.path(
    #       dir.figures, 
    #       paste(sep = "", "DESeq2_vs_edgeR_Venn",
    #             "_", s, "_", thresholds[s], 
    #             "_", prefix["comparison"], ".pdf")))
    vennDiagram(venn.counts.one.threshold, cex=1, main=paste(sep=" ", cond1, "vs", cond2, "; ", s, "threshold = ", thresholds[s]))
   silence <- dev.off()
  }
  
  
  ################################################################
  ## Compare counts per million reads (CPMs) between samples. 
  ## This is just to get an intuitive idea, since CPMs are 
  ## not recommended for diffferential detection.
  verbose("\t\tmean CPM plot", 2)
  pdf(file=paste(sep = "", prefix["comparison_figure"], "CPM_plot",  ".pdf"))
  plot(result.table[,c("cpm1.mean", "cpm2.mean")], 
       log="xy",
       main = "Mean counts per million reads (log scales)",
       xlab=paste(cond1),
       ylab=paste(cond2),
       col=gene.colors,
       panel.first=grid(lty="solid", col="#DDDDDD"))
  ## Plot genes on the top layer to highlight them
  points(result.table[both,c("cpm1.mean", "cpm2.mean")], col=gene.palette["both"])
  points(result.table[DESeq2.only,c("cpm1.mean", "cpm2.mean")], col=gene.palette["DESeq2.only"])
  points(result.table[edgeR.only,c("cpm1.mean", "cpm2.mean")], col=gene.palette["edgeR.only"])
  abline(a=0, b=1)
  legend("topleft", 
         pch=1, cex=1, 
         bg="white", bty="o",
         col=gene.palette[c("both", "DESeq2.only", "edgeR.only", "none")], 
         legend = c(paste(sum(both, na.rm=TRUE), "both"),
                    paste(sum(DESeq2.only, na.rm=TRUE), "DESeq2 only"),
                    paste(sum(edgeR.only, na.rm=TRUE), "edgeR only"),
                    paste(sum(none, na.rm=TRUE), "none")))
  silence <- dev.off()
  
  
  ################################################################
  ## Draw MA plot with CPMs
  verbose("\t\tCPM MA plot", 2)
  pdf(file=paste(sep = "", prefix["comparison_figure"], "_CPM_MA_plot.pdf"))
  #pdf(file=file.path(dir.figures, paste(sep = "", "CPM_MA_plot_", prefix["comparison"], ".pdf")))
  plot(result.table[,c("A", "M")],
       main = paste(sep=" ", cond1, "vs", cond2, ": CPMs MA plot"),
       xlab=paste(sep="", "A = log2(", cond1, "*", cond2, ")/2"),
       ylab=paste(sep="", "M = log2(", cond1, "/", cond2, ")"),
       col=gene.colors,
       panel.first=grid(lty="solid", col="#DDDDDD"))
  ## Plot genes on the top layer to highlight them
  points(result.table[both,c("A", "M")], col=gene.palette["both"])
  points(result.table[DESeq2.only,c("A", "M")], col=gene.palette["DESeq2.only"])
  points(result.table[edgeR.only,c("A", "M")], col=gene.palette["edgeR.only"])
  abline(h=0)
  legend("bottomright", 
         pch=1, cex=1, 
         bg="white", bty="o",
         col=gene.palette[c("both", "DESeq2.only", "edgeR.only", "none")], 
         legend = c(paste(sum(both, na.rm=TRUE), "both"),
                    paste(sum(DESeq2.only, na.rm=TRUE), "DESeq2 only"),
                    paste(sum(edgeR.only, na.rm=TRUE), "edgeR only"),
                    paste(sum(none, na.rm=TRUE), "none")))
  silence <- dev.off()
  
 
  ## Volcano plot
#   for (deg.tool in deg.tools) {
#     criterion <- "padj"
#     for (criterion in c("padj", "evalue")) {
#       verbose(paste(sep=" ", "\t\t", deg.tool, criterion, "Volcano plot"), 3)
#       pdf(file=paste(sep="_", prefix[paste(sep="", deg.tool, "_figure")], 
#                      criterion, thresholds[criterion], "VolcanoPlot.pdf"))
#       control.col <- paste(sep=".", deg.tool, criterion)
#       effect.col <- paste(sep=".", deg.tool, "log2FC")
#       VolcanoPlot(na.omit(result.table[,c(control.col, effect.col)]), 
#                   control.type=control.col, alpha = thresholds[criterion], 
#                   effect.size.col=effect.col, xlab="log2(fold-change)", effect.threshold=log2(thresholds["FC"]),
#                   main=paste(sep=" ", cond1, "vs", cond2, "; ", deg.tool, " ", criterion, "volcano plot"),
#                   sort.by.pval = TRUE, plot.points = TRUE,
#                   legend.corner = "topleft")
#       quiet <- dev.off()
#     }
#    }  
  
  ################################################################
  ## Store the Differentally expressed genes in a two-column file
  # names(result.table)
  DEG.columns  <- c("DEG_edgeR_and_DESeq2", "edgeR.DEG", "DESeq2.DEG", "padj_edgeR_and_DESeq2")
  # column <- DEG.columns[3]
  deg.id.lists <- list()
  column <- DEG.columns[1]
  for (column in DEG.columns) {
    is.DEG <- !is.na(result.table[,column]) & result.table[,column] == 1
    # table(is.DEG)
    DEG.ids <- na.omit(row.names(result.table[is.DEG,]))
    deg.id.lists[[column]] <- DEG.ids
    # length(DEG.ids)
    # table(is.na(DEG.ids))
    DEG.genes <- data.frame(gene.id=DEG.ids, "cluster"=rep(column, times=length(DEG.ids)))
    if (column == DEG.columns[1]) {
      cluster.table <- DEG.genes
    } else {
      cluster.table <- rbind(cluster.table, DEG.genes)
    }
  }
  # View(cluster.table)
  # table(cluster.table$cluster)
  # grep (cluster.table$gene.id, pattern="NA", value=TRUE)
  ## Save gene clusters
  cluster.file <- paste(sep = "", 
                        prefix["comparison_file"], 
                        "_", suffix.deg, "_DESeq2_edgeR_clusters.tab")
  write.table(x = cluster.table, row.names = FALSE,
              file = cluster.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tCluster file\t", cluster.file), 1)
  
  ################################################################
  ## Functional enrichment analysis
  run.functional.enrichment <- FALSE
  if (run.functional.enrichment & exists("org.db") & !is.na(org.db) & !is.null(org.db) & exists("gene.info.rsat")) {
    library(org.db, character.only = TRUE)
    verbose("Starting the analysis of functional enrichment")
    
    ## Convert IDs to entrez IDs
    verbose("Converting gene names to Entrez IDs")
    all.gene.ids <- row.names(result.table)
    if (sum(!is.na(gene.info$entrez.id)) == 0) {
      gg <- bitr(row.names(result.table), fromType="SYMBOL", toType="ENTREZID", annoDb=org.db, drop=FALSE)
      row.names(gg) <- gg$SYMBOL
      # dim(gg)
      # dim(result.table)
      # head(gg)
      all.entrez.ids <- gg[all.gene.ids, "ENTREZID"]
    } else {
      all.entrez.ids <- as.vector(gene.info[all.gene.ids,"entrez.id"])      
    }
    names(all.entrez.ids) <- gene.info[all.gene.ids,"gene_id"]
    # table(!is.na(all.entrez.ids))
    all.entrez.ids <- unique(na.omit(all.entrez.ids))
    # length(all.entrez.ids)
    
    ## Select the selection columns on wchich enrichment analysis will be performed.
    ## We select each single-score threshold column, as well as the combined thresholds ("DEG" columns).
    # names(result.table)
    geneset.selection.columns <- vector()
    for (deg.tool in deg.tools) {
      for (s in names(thresholds)) {
        geneset.selection.columns <- append(geneset.selection.columns, paste(sep="", deg.tool, ".", s, "_", thresholds[s]))
      }      
      geneset.selection.columns <- append(geneset.selection.columns, paste(sep="", deg.tool, ".DEG"))
    }
    # print(geneset.selection.columns)
    
    col <- "edgeR.DEG"
    for (col in geneset.selection.columns) {
      verbose(paste("Gene selection column", col))
      
      geneset.ids <- as.vector(as.matrix(all.gene.ids[result.table[,col] == 1]))
      geneset <- all.entrez.ids[geneset.ids]
      # table(is.na(geneset))
      geneset <- geneset[!is.na(geneset)]
      gene.nb <- length(geneset)
      verbose(paste("\t", gene.nb, "selected genes"))
      if (gene.nb == 0) {
        verbose(paste("\t", "Skipping: not a single gene selected"))
        next
      }
      
      # length(geneset)
      verbose(paste("\tGeneFunctional enrichment"))
      enrich.result <- functional.enrichment(geneset=geneset,
                                             allgenes=all.entrez.ids,
                                             db=org.db,
                                             ontology="BP",
                                             thresholds = c("evalue"=1, "qvalue"=0.05),
                                             select.positives=FALSE,
                                             run.GOstats = TRUE,
                                             run.clusterProfiler = FALSE,
                                             organism.names=organism.names,
                                             plot.adjust = TRUE)
      
      ## Select GO Biological process result table
      go.bp.table <- enrich.result$go.bp.table
      go.bp.table.positive <- go.bp.table[go.bp.table$positive==1,]
      
      # View(enrich.result$go.bp.table)
      ## Save full result table
      go.bp.file <- paste(sep = "", prefix["comparison_file"], 
                          "_", suffix.deg,
                          "_", col,
                          "_GOstats_all.tab")
      write.table(x = go.bp.table, row.names = FALSE,
                  file = go.bp.file, sep = "\t", quote=FALSE)
      verbose(paste(sep="", "\tGOstats BP over-representation\t", 
                    nrow(go.bp.table), " rows\t", go.bp.file), 1)
      
      ## Save positive associations in a separate file
      go.bp.file.positive <- paste(sep = "", prefix["comparison_file"], 
                                   "_", suffix.deg,
                                   "_", col,
                                   "_GOstats_positive.tab")
      write.table(x = go.bp.table.positive, row.names = FALSE,
                  file = go.bp.file.positive, sep = "\t", quote=FALSE)
      verbose(paste(sep="", "\tGOstats BP over-representation\t", 
                    nrow(go.bp.table.positive), " rows\t", go.bp.file.positive), 1)
      
      #       kk <- enrichKEGG(gene = geneset,
      #                        organism = "eco",
      #                        pvalueCutoff = 1,
      #                        readable = TRUE)
      #       head(summary(kk))
      #       
      ## Build a custom GO map from 2-column data.frame with GO, gene (entrez ID)
      if ((exists("go.map.file")) & (exists("go.description.file"))) {
        ## Read description of GO terms
        go.description <- read.delim(go.description.file, header=TRUE, row.names=1)
        # dim(go.description)
        # names(go.description)
        # head(go.description)
        
        ## Read the GO-gene annotation table (tab-delimited file)
        gomap.frame <- read.delim(go.map.file, header=FALSE)
        names(gomap.frame) <- c("gene.id", 
                                "GO.id", 
                                "gene.name", 
                                "transcript.name", 
                                "protein.id", 
                                "entrez.id")
        gomap.frame$GO.descr <- go.description[as.vector(gomap.frame$GO.id), "GO.Term"] ## Load GO term descriptions in GO-gene table
        # dim(gomap.frame)
        # head(gomap.frame)
        
        #go.annot <- buildGOmap(gomap.frame[,c(2,6)])
        go.enricher.res = enricher(geneset, 
                                   TERM2GENE=gomap.frame[, c("GO.id", "entrez.id")], 
                                   TERM2NAME=gomap.frame[, c("GO.id", "GO.descr")],
                                   pAdjustMethod = "BH", 
                                   minGSSize = 0,
                                   pvalueCutoff = 1,
                                   qvalueCutoff=1)
        go.enricher.table <- data.frame(go.enricher.res@result)
        go.enricher.table <- complete.enrich.table(go.enricher.table, pvalue.column = "pvalue")
        # dim(go.enricher.table)
        # View(go.enricher.table)
        ## Select positive rows
        go.enricher.table.positive <- go.enricher.table[go.enricher.table$positive==1,]
        
        ## Save result table
        go.enricher.file <- paste(sep = "", prefix["comparison_file"], 
                                  "_", suffix.deg,
                                  "_", col,
                                  "_GO_enricher_all.tab")
        write.table(x = go.enricher.table, row.names = FALSE,
                    file = go.enricher.file, sep = "\t", quote=FALSE)
        verbose(paste(sep="", "\tGO clusterProfiles::enricher\t", 
                      nrow(go.enricher.table), " rows\t", go.enricher.file), 1)
        
        ## Save result table
        go.enricher.file.positive <- paste(sep = "", prefix["comparison_file"], 
                                           "_", suffix.deg,
                                           "_", col,
                                           "_GO_enricher_positive.tab")
        write.table(x = go.enricher.table.positive, row.names = FALSE,
                    file = go.enricher.file.positive, sep = "\t", quote=FALSE)
        verbose(paste(sep="", "\tGO clusterProfiles::enricher\t", 
                      nrow(go.enricher.table.positive), " rows\t", go.enricher.file.positive), 1)
      }
    }
    
    
    #   library("GenomicFeatures")
    #   gtf.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.gtf"
    #   txdb <- makeTxDbFromGFF(file=gtf.file,
    #                           organism="Escherichia coli",
    #                           #                         genome="Escherichia_coli_str_k_12_substr_mg1655",
    #                           dataSource="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-28/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/")
    #   seqlevels(txdb)
    #   genome(txdb) <- "Escherichia_coli_str_k_12_substr_mg1655"
    #   
    #   # Build specific GO map using GFF file
    #   #   library(biomaRt)
    #   #   Mtb <- useMart(biomart="fungi_mart_28",
    #   #                  dataset="scerevisiae")
    #   #   listDatasets(mart = Mtb)
    #   #   
    #   #   mmart <-
    #   #     makeTxDbFromBiomart(biomart = "fungi_mart_28",
    #   #                         dataset = "scerevisiae_gene_ensembl")
    #   #   listDatasets(biomart = "fungi_mart_28")
    #   #   library("biomaRt")
    #   #   listMarts()[1] ## The bacterial Mart has disappeared in version 28 !!!!
    #   library("clusterProfiler")
    #   library("org.EcK12.eg.db")
    #   eg <- bitr(gene.ids, fromType="SYMBOL", toType="ENTREZID", annoDb="org.EcK12.eg.db")
    #   
  }
  
  ################################################################
  ## Summarise results of the current analysis
  
  ## Instantiate a data frame for the current analysis
  current.summary <- data.frame(
    "analysis"=paste(sep="", prefix["comparison"]),
    "cond1" = cond1, "cond2" = cond2)
  
  ## DESeq2 results
  for (s in names(thresholds)) {
    for (deg.tool in deg.tools) {
      selection.column <- paste(sep="", deg.tool, ".", s, "_", thresholds[s])
      current.summary[, selection.column] <- sum(result.table[, selection.column], na.rm=TRUE)
    }
    
    for (combination in c("edgeR_and_DESeq2",
                          "edgeR_not_DESeq2",
                          "DESeq2_not_edgeR",
                          "none")) {
      
      selection.column <- paste(sep="_", s, combination)
      current.summary[, selection.column] <- sum(result.table[, selection.column], na.rm=TRUE)
    }
    
    ## Export gene lists
    
  }
  
  if (i == 1) {
    summary.per.analysis <- current.summary
  } else {
    summary.per.analysis <- rbind(summary.per.analysis, current.summary)
  } 
}

########################################################################
## Export summary table
verbose("Exporting summary table", 1)
summary.file <- paste(sep="", prefix["general.file"], "_summary_per_analysis.tab")
write.table(x = t(summary.per.analysis), row.names = TRUE, col.names=FALSE,
            file = summary.file, sep = "\t", quote=FALSE)
verbose(paste(sep="", "\tSummary per analysis\t", summary.file), 1)

########################################################################
## Save an image of memory in order to reload the whole analysis without 
## re-computing everything. 
## This memory image can also be directly loaded into memory on another computer.
if (save.image) {
  image.file <- paste(sep="", prefix["general.file"], "_memory_image.RData")
  verbose(paste("Working directory", getwd()), 1)
  verbose(paste("Saving memory image", image.file), 1)
  save.image(file=image.file, compress=TRUE)
} else {
  verbose("Skipping memory image saving")
}


