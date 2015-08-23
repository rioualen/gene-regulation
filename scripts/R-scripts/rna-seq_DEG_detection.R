################################################################
## R code to detect differentially expressed genes with a combination 
## of several  BioConductor libraries, including edgeR and DESeq2.
##
## Authors: Jeanne Chèneby & Justine Long (Initial script), revised and completed by Jacques van Helden
## First version: 2015-07
## Last update: 2015-08-23
################################################################

library("edgeR", warn.conflicts = FALSE, quietly=TRUE)
library("DESeq2", warn.conflicts = FALSE, quietly=TRUE, verbose=FALSE)
library("limma", warn.conflicts = FALSE, quietly=TRUE) ## Required for vennCounts and vennDiagram
library(gplots, warn.conflicts = FALSE, quietly=TRUE) ## Required for heatmaps.2
library(RColorBrewer, warn.conflicts = FALSE, quietly=TRUE)
library("GenomicFeatures")

## Load library of functions for differential analysis
dir.fg <- "~/fg-chip-seq/"
source(file.path(dir.fg, "scripts/R-scripts/deg_lib.R"))


verbosity <- 1
deg.tools <- c("edgeR", "DESeq2")

## Elements that should be added to the parameters
org.db <- "org.EcK12.eg.db" ## Should be added to parameters
gene.info.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2_gene_info.tab"
organism.name <- "Escherichia coli"
gtf.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.gtf"
gtf.source <- "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-28/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/"
#   pet.gene <- "b2531"

## The only argument is the file containing all the parameters for the analysis
if (is.na(commandArgs(trailingOnly = FALSE)[6])) {
  #   stop("Parameter file must be passed as command-line argument.")
  ## TEMPORARY FOR DEBUGGING: 
  setwd("~/BeatriceRoche/")
  r.params.path <- "results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_params.R"  
} else {
  r.params.path <- commandArgs(trailingOnly = FALSE)[6]  
}
verbose(paste("Reading parameters", r.params.path), 1)
source(r.params.path)
setwd(dir.main)

################################################################
## Check parameters + define additional ones

## Choose default value for thresholds if they were not defined in the config file
if (!exists("thresholds")) {
  thresholds <- c("padj"=0.05, "evalue"=1, "FC"=1.5)
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

## Define a specific color for each distinct condition
exp.conditions <- unique(sample.conditions) ## Set of distinct conditions
cols.conditions <- brewer.pal(length(exp.conditions),"Dark2")
names(cols.conditions) <- exp.conditions

## Define a color per sample according to its condition
cols.samples <- cols.conditions[sample.conditions]
names(cols.samples) <- sample.ids

## Read the design file, which indicates the anlayses to be done.
## Each row specifies one differential expression analysis, which 
## consists in comparing two conditions. 
verbose("Reading design", 1)
design <- read.delim(design.file, sep="\t", 
                     comment=";", header = TRUE, row.names=NULL)

## Prefix for output files concerning the whole count table (all samples together)
## prefix["general.file"] <- sub(pattern = ".tab", replacement="", all.counts.table)
prefix <- vector()
prefix["general.file"] <- file.path(dir.DEG, suffix.deg)

## Read the count table
verbose("Loading count table", 1)
all.counts <- read.delim(all.counts.table, row.names=1, sep="\t")
# names(all.counts)

verbose("Computing count-derived metrics (log-transformed)", 1)

## Statistics on reads that could not be mapped to genes for different reasons (intergenic, ambiguous, not unique, ...)
all.counts.not.mapped <- all.counts[grep(pattern = "^__", x = row.names(all.counts)), ]
# dim(all.counts.not.mapped)
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
  #   library(rtracklayer)  # for the import() function
  #   gr <- import(gtf.file)
  txdb <- makeTxDbFromGFF(file=gtf.file,
                          organism=organism.name,
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
## TEMPORARY: load gene information from a tab-delimited file generated with RSAT, 
## because for bacterial genomes there is no obvious way to obtain EntrezIDs.
if (exists("gene.info.file")) {
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
  # View(gene.info)
} 

## Export the gene information table to keep a trace of what has been used. 
verbose("Exporting gene information table", 1)
gene.info.out <- paste(sep="", prefix["general.file"], "_gene_descriptions.tab")
write.table(x = gene.info, row.names = FALSE,
            file = gene.info.out, sep = "\t", quote=FALSE)
verbose(paste(sep="", "\tGene info table\t", gene.info.out), 1)


################################################################
## Compute sample-wise statistics on mapped counts
################################################################
verbose("Computing statistics per sample", 1)

stats.per.sample <- data.frame(
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

## Count number and the fraction of samples with counts below the mean. 
## This shows the impact of very large counts: in my test samples, 
## 85% of the samples have a value below the mean (i.e. the mean is at the percentile 85 !)
stats.per.sample$below.mean <- apply(t(all.counts.mapped) < stats.per.sample$mean, 1, sum, na.rm=TRUE)
stats.per.sample$fract.below.mean <- stats.per.sample$below.mean/nrow(all.counts.mapped)

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
cpms.median <- cpm(all.counts.mapped.epsilon, lib.size = stats.per.sample$median)    ## Counts per million reads, normalised by sample-wise median count
cpms <- cpms.perc75 ## Choose one normalization factor for the CPMs used below
cpms.log10 <- log10(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts
cpms.log2 <- log2(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts

################################################################
## Draw some generic plots
################################################################

verbose("Drawing generic plots from the whole count table", 1)


## Plot the impact of the normalization factor (library sum , median or percentile 75)
png(file= file.path(dir.DEG, paste(sep = "", "CPM_libsum_vs_median_vs_perc75.png")), 
    width=1000, height=1000)
par.ori <- par() ## Save original plot parameters
cols.counts <- as.data.frame(matrix(cols.samples, nrow=nrow(all.counts.mapped), ncol=ncol(all.counts.mapped), byrow = TRUE))
colnames(cols.counts) <- names(all.counts.mapped)
rownames(cols.counts) <- rownames(all.counts.mapped)
plot(data.frame("libsum" = as.vector(as.matrix(cpms.libsum)),
                "median" = as.vector(as.matrix(cpms.median)),
                "perc75" = as.vector(as.matrix(cpms.perc75))),
     col=as.vector(as.matrix(cols.counts)))
quiet <- dev.off()

## Plot some sample-wise statistics
pdf(file= file.path(dir.DEG, paste(sep = "", "sample_statistics_plots.pdf")), width=10, height=10)
par(mar=c(5,5,1,1)) ## adpt axes
par(mfrow=c(2,2))
## Median versus mean
plot(stats.per.sample[,c("mean", "median")], 
     panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
     las=1, col=cols.samples)

## First versus third quartile
plot(stats.per.sample[,c("perc25", "perc75")], 
     panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
     las=1, col=cols.samples)

## Sum versus third quartile. 
plot(stats.per.sample[,c("sum", "perc75")], 
     panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
     las=1, col=cols.samples)

## Mean versus third quartile. 
plot(stats.per.sample[,c("mean", "perc75")], 
     panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
     las=1, col=cols.samples)
par(mfrow=c(1,1))
quiet <- dev.off()

################################################################
## Boxplots of raw counts and derived counting measures
par(mar=c(5,7,4,1)) ## adapt axes

## Boxplot of raw counts
pdf(file= file.path(dir.DEG, paste(sep = "", "sample_boxplots_raw_counts.pdf")), width=7, height=7)
boxplot(all.counts.mapped, horizontal=TRUE, col=cols.samples,
        xlab="Raw counts",
        main="Box plots per sample: raw counts", las=1)
quiet <- dev.off()

## Boxplot of log10-transformed counts
pdf(file= file.path(dir.DEG, paste(sep = "", "sample_boxplots_log10_counts.pdf")), width=7, height=7)
boxplot(all.counts.mapped.log10, horizontal=TRUE, col=cols.samples,
        xlab="log10(counts)",
        main="Box plots per sample: log10(counts)", las=1)
quiet <- dev.off()

## Boxplot of CPMs
pdf(file= file.path(dir.DEG, paste(sep = "", "sample_boxplots_CPM.pdf")), width=7, height=7)
boxplot(cpms, horizontal=TRUE, col=cols.samples,
        xlab="CPM",
        main="Box plots per sample: counts per million reads (CPM)", las=1)
quiet <- dev.off()

## Boxplot of log10-transformed CPMs
pdf(file= file.path(dir.DEG, paste(sep = "", "sample_boxplots_log10_CPM.pdf")), width=7, height=7)
boxplot(cpms.log10, horizontal=TRUE, col=cols.samples,
        xlab="log10(CPM)",
        main="Box plots per sample: counts per million reads (CPM)", las=1)
quiet <- dev.off()

par <- par.ori ## Restore original plot parameters


## Draw sample correlation heatmaps for the raw read counts
pdf(file=paste(sep="", prefix["general.file"],"_sample_correl_heatmap_counts.pdf"))
hm <- heatmap.2(as.matrix(cor(all.counts.mapped)),  scale="none", trace="none", 
                main="Correlation between raw counts",
                col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
quiet <- dev.off()

## Draw sample correlation heatmaps for CPM
#heatmap(cor(all.counts.mapped), scale = "none")
pdf(file=paste(sep="", prefix["general.file"],"_sample_correl_heatmap_cpms.pdf"))
hm <- heatmap.2(as.matrix(cor(cpms)),  scale="none", trace="none", 
                main="Correlation between CPM",
                col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
quiet <- dev.off()

## Plot the first versus second components of samples
cpms.pc <- prcomp(t(cpms))
pdf(file=paste(sep="", prefix["general.file"],"_CPM_PC1-PC2.pdf"))
plot(cpms.pc$x[,1:2], panel.first=grid(), type="n", main="First components from PCA-transformed CPMs")
text(cpms.pc$x[,1:2], labels = sample.conditions, col=cols.samples)
quiet <- dev.off()


################################################################
## Analyse between-replicate reproducibility
################################################################
verbose("Plotting betwen-replicate comparisons", 1)
cond <- "cyay"
for (cond in exp.conditions) {
  verbose(paste(sep="", "\tcondition\t", cond), 2)
  
  max.rep.to.plot <- 5 ## Restrict the number of replicates to plot, for the sale of readability
  
  ## Create a specific result directory for this condition
  dir.condition <- file.path(dir.DEG, "per_condition", cond)
  dir.create(path = dir.condition, showWarnings = FALSE, recursive = TRUE)
  
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
  
  
  ## Plot pairwise comparisons between replicates
  pdf(file= file.path(dir.condition, paste(sep = "", "between-replicate_counts_plot_", cond, ".pdf")), width=10, height=10)
  plot(current.counts[,1:min(max.rep.to.plot, nrep)], log="xy", col=cols.conditions[cond], 
       main=paste(cond, " ; raw counts per replicate (log scale)"))
  dev.off()
  
  ## Plot mean versus variance of raw counts for the current condition
  pdf(file= file.path(dir.condition, paste(sep = "", "counts_variance-mean_plot_", cond, ".pdf")), width=10, height=10)
  plot(current.counts.mean, 
       current.counts.var, 
       log="xy", col=cols.conditions[cond], 
       panel.first=grid(lty="solid", col="#DDDDDD"),
       xlab=paste("mean counts for condition", cond),
       ylab=paste("variance of counts for condition", cond),
       main=paste(cond, " ; Counts variance/Mean plot"))
  abline(a=0,b=1, lty="dashed", col="green", lwd=2) ## Milestone for Poisson distributions: var = mean
  quiet <- dev.off()
  
  ## Plot mean versus variance of CPMs for the current condition
  pdf(file= file.path(dir.condition, paste(sep = "", "CPM_variance-mean_plot_", cond, ".pdf")), width=10, height=10)
  plot(current.cpm.mean, 
       current.cpm.var, 
       log="xy", col=cols.conditions[cond], 
       panel.first=grid(lty="solid", col="#DDDDDD"),
       xlab=paste("CPM mean for condition", cond),
       ylab=paste("CPM variance for condition", cond),
       main=paste(cond, " ; CPM variance/Mean plot"))
  abline(a=0,b=1, lty="dashed", col="green", lwd=2) ## Milestone for Poisson distributions: var = mean
  quiet <- dev.off()
  
}

################################################################
## Run differential expression analysis
################################################################
verbose("Starting differential analysis", 1)

## Iterate over analyses
i <- 1
for (i in 1:nrow(design)) {
  
  ## Identify samples for the first condition
  cond1 <- as.vector(design[i,1])  ## First condition for the current comparison
  samples1 <- sample.ids[sample.conditions == cond1]
  
  ## Identify samples for the second condition
  cond2 <- as.vector(design[i,2])  ## Second condition for the current comparison
  samples2 <- sample.ids[sample.conditions == cond2]
  
  verbose(paste(sep="", "\tDifferential analysis\t", i , "/", nrow(design), "\t", cond1, " vs ", cond2), 1)
  
  ## Create a specific result directory for this differential analysis
  prefix["comparison"] <- paste(sep="_", cond1, "vs", cond2)
  dir.analysis <- file.path(dir.DEG, paste(sep="", prefix["comparison"]))
  dir.create(path = dir.analysis, showWarnings = FALSE, recursive = TRUE)
  dir.figures <- file.path(dir.analysis, "figures")
  dir.create(path = dir.figures, showWarnings = FALSE, recursive = TRUE)
  prefix["comparison_file"] <- file.path(dir.analysis, prefix["comparison"])
  prefix["comparison_figure"] <- file.path(
    dir.figures, 
    paste(sep="", prefix["comparison"], "_",  suffix.deg))
  
  ## Select counts for the samples belonging to the two conditions
  current.samples <- c(samples1, samples2)
  current.counts <- all.counts.mapped[,current.samples]
  # dim(current.counts)  ## For test
  
  if (sum(!names(current.counts) %in% sample.ids) > 0) {
    stop("Count table contains column names without ID in sample description file.")
  }
  
  ## Define conditions and labels for the samples of the current analysis
  current.conditions <- sample.conditions[current.samples]
  current.labels <- paste(current.conditions, names(current.counts), sep="_")
  
  ## Initiate a result table with the CPMs and derived statistics
  all.gene.ids <- row.names(cpms)
  result.table <- data.frame("gene_id" = all.gene.ids,
                             "name"=gene.info[all.gene.ids,"name"])
  row.names(result.table) <- all.gene.ids
  result.table$entrez.id <- gene.info[all.gene.ids,"entrez.id"]
  result.table$description <- gene.info[all.gene.ids,"description"]
  #   result.table$mean.cpm1 <- apply(cpms[,samples1],1, mean)
  #   result.table$mean.cpm2 <- apply(cpms[,samples2],1, mean)
  #   result.table$A <- log2(result.table$mean.cpm1*result.table$mean.cpm2)/2
  #   result.table$M <- log2(result.table$mean.cpm1/result.table$mean.cpm2)
  # dim(result.table)
  
  ## Tag genes detected in less than min.rep samples, which is defined as 
  ## the minimal number of replicates per condition.
  min.rep <- min(length(samples1), length(samples2))
  result.table$undetected <- rowSums(current.counts > 1) < min.rep
  #current.counts <- current.counts[rowSums(current.counts > 1) >= min.rep,]
  # dim(current.counts)
  
  result.table$cpm.mean <- apply(cpms,1, mean)
  result.table$cpm1.mean <- apply(cpms[,samples1],1, mean)
  result.table$cpm2.mean <- apply(cpms[,samples2],1, mean)
  result.table$A = log2(result.table$cpm1.mean*result.table$cpm2.mean)/2
  result.table$M = log2(result.table$cpm1.mean/result.table$cpm2.mean)
  result.table$cpm.median <- apply(cpms,1, median)
  result.table$cpm1.median <- apply(cpms[,samples1],1, median)
  result.table$cpm2.median <- apply(cpms[,samples2],1, median)
  result.table$cpm.min <-  apply(cpms,1, min)
  result.table$cpm1.min <- apply(cpms[,samples1],1, min)
  result.table$cpm2.min <- apply(cpms[,samples2],1, min)
  result.table$cpm.max <-  apply(cpms,1, max)
  result.table$cpm1.max <- apply(cpms[,samples1],1, max)
  result.table$cpm2.max <- apply(cpms[,samples2],1, max)
  result.table$cpm.sd <-  apply(cpms,1, sd)
  result.table$cpm1.sd <- apply(cpms[,samples1],1, sd)
  result.table$cpm2.sd <- apply(cpms[,samples2],1, sd)
  result.table$cpm.var <-  apply(cpms,1, var)
  result.table$cpm1.var <- apply(cpms[,samples1],1, sd)
  result.table$cpm2.var <- apply(cpms[,samples2],1, sd)
  # View(result.table)
  
  ################################################################
  ## DESeq2 analysis
  ################################################################
  
  verbose("\t\tDESeq2 analysis", 2)
  
  ## Path prefix to save DESeq2 result files
  prefix["DESeq2_file"] <- paste(sep="", prefix["comparison_file"], "_", suffix.DESeq2)
  prefix["DESeq2_figure"] <- paste(sep="", prefix["comparison_figure"], "_", suffix.DESeq2)
  
  ## Create a DESeqDataSet object from the count table + conditions
  condition <- as.factor(as.vector(current.conditions))
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
    thresholds=c("padj"=0.05, "evalue"=1, "FC"=1.5),
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
  write.table(x = edger.result.table, 
              row.names = FALSE, file = edger.result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tedgeR result file\t", edger.result.file), 1)
  
  ################################################################
  ## Export edgeR plots
  
  verbose("\t\tGenerating edgeR plots", 2)
  
  ## Smear plot
  verbose("\t\t\tedgeR Smear plot", 3)
  edger.deg <- rownames(edger.tt$table)[edger.tt$table$FDR < thresholds["padj"]] ## List of differentially expressed genes reported by edgeR
  pdf(file=paste(sep="", prefix["DESeq2_figure"], "_plotSmear.pdf"))
  plotSmear(d, de.tags = edger.deg)
  quiet <- dev.off()
  
  ## MDS plot (Multidimensional scaling plot of distances between gene expression profiles)
  verbose("\t\t\tedgeR MDS plot", 3)
  pdf(file=paste(sep="", prefix["DESeq2_figure"], "_plotMDS.pdf"))
  #  pdf(file=file.path(dir.figures, paste(sep="", "edgeR_MDS_plot_", prefix["comparison"], ".pdf")))
  plotMDS(d, labels=current.labels, 
          col=c("darkgreen","blue")[factor(sample.conditions[names(current.counts)])])
  quiet <- dev.off()
  
  # Mean-variance relationship
  verbose("\t\t\tedgeR mean-variance relationship", 3)
  pdf(file=paste(sep="", prefix["DESeq2_figure"], "_plotMeanVar.pdf"))
  plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
  quiet <- dev.off()
  
  
  # BCV (Biological Coefficient of Variation)
  verbose("\t\t\tedgeR BCV plot", 3)
  pdf(file=paste(sep="", prefix["DESeq2_figure"], "_plotBCV.pdf"))
  #pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotBCV_", prefix["comparison"], ".pdf")))
  plotBCV(d)
  quiet <- dev.off()
  
  ################################################################
  ## Compare gene selections by edgeR and DESeq2
  verbose("\t\tComparing DESeq2 and edgeR results", 2)    
  ## Add tags indicating the combinations of edgeR/DESeq2 selections
  for (s in c("padj", "evalue", "FC")) {
    selection.column <- paste(sep="_", s, thresholds[s])
    
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
  write.table(x = result.table, row.names = FALSE,
              file = result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tMerged result file\t", result.file), 1)
  
  ## Define point colors according to the test results
  both <- result.table$padj_edgeR_and_DESeq2 == 1
  edgeR.only <- result.table$padj_edgeR_not_DESeq2 == 1
  DESeq2.only <- result.table$padj_DESeq2_not_edgeR == 1
  none <- result.table$padj_none
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
  dev.off()
  
  
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
    dev.off()
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
  dev.off()
  
  
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
  dev.off()
  
  
  ################################################################
  ## Functional enrichment analysis
  if (exists("org.db") & !is.na(org.db) & !is.null(org.db) & exists("gene.info.rsat")) {
    
    all.gene.ids <- row.names(result.table)
    
    ## Convert IDs to entrez IDs using custom annotation table
    entrez.ids <- as.vector(gene.info[all.gene.ids,"entrez.id"])
    names(entrez.ids) <- gene.info[all.gene.ids,"gene_id"]
    
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
      geneset <- all.gene.ids[result.table[,col] == 1]
      go.bp.table <- gostat.overrepresentation(geneset=na.omit(entrez.ids[geneset]), 
                                               allgenes=entrez.ids, 
                                               db=org.db, evalue.filter=TRUE,
                                               verbosity=2)
      # View(go.bp.table)
      ## Save result table
      go.bp.file <- paste(sep = "", prefix["comparison_figure"], 
                          "_", suffix.deg,
                          "_", col,
                          "_GOstats.tab")
      write.table(x = go.bp.table, row.names = FALSE,
                  file = go.bp.file, sep = "\t", quote=FALSE)
      verbose(paste(sep="", "\t\tGO over-representation analysis file\t", go.bp.file), 1)
      
      
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

## Export summary table
verbose("Exporting summary table", 1)
summary.file <- paste(sep="", prefix["general.file"], "_summary_per_analysis.tab")
write.table(x = summary.per.analysis, row.names = FALSE,
            file = summary.file, sep = "\t", quote=FALSE)
verbose(paste(sep="", "\tSummary per analysis\t", summary.file), 1)


