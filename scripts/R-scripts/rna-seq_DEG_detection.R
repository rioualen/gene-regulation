################################################################
## R code to detect differentially expressed genes with a combination 
## of several  BioConductor libraries, including edgeR and DESeq2.
##
## Authors: Jeanne Chèneby & Justine Long & Jacques van Helden
## First version: 2015-07
################################################################

library("edgeR", warn.conflicts = FALSE, quietly=TRUE)
library("DESeq2", warn.conflicts = FALSE, quietly=TRUE, verbose=FALSE)
library("limma", warn.conflicts = FALSE, quietly=TRUE) ## Required for vennCounts and vennDiagram
library(gplots, warn.conflicts = FALSE, quietly=TRUE) ## Required for heatmaps.2
library(RColorBrewer, warn.conflicts = FALSE, quietly=TRUE)

verbosity <- 1

#' @title Display messages at a given verbosity level
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Display messages depending on user-defined verbosity level
#'
#' @details
#' First version: 2015-03. 
#' Last modification: 2015-03. 
#'
#' @param verbosity   Level of verbosity above which the message should be printed.
#' @param print.date=TRUE   Print date and time
#'
#' @examples
#'
#' verbosity <- 1 ## Define level of verbosity
#'
#' ## This message will be printed because the level is <= verbosity
#' verbose("This is printed", 1)
#'
#' ## This message will not be printed because the verbosity is inferior to the specified level
#' verbose("This is not printed", 2)
#'
#' @export
verbose <- function(message.content,
                    level=1,
                    print.date=TRUE) {
  if (!exists("verbosity")) {
    verbosity <- 1
  }
  if (verbosity >= level) {
    if (print.date) {
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t", message.content)
    } else {
      message(message.content)
    }
  }
}

verbose("Reading parameters", 1)

## Define a color palette for heatmaps. I like this Red-Blue palette because 
## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
## - it can be seen by people suffering from red–green color blindness.
cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))

## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
epsilon <- 0.01

## The only argument is the file containing all the parameters for the analysis
r.params.path <- commandArgs(trailingOnly = FALSE)[6]
## TEMPORARY FOR DEBUGGING: 
## setwd("~/BeatriceRoche/")
#  r.params.path <- "results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_params.R"
source(r.params.path)
setwd(dir.main)

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
## all.prefix <- sub(pattern = ".tab", replacement="", all.counts.table)
all.prefix <- file.path(dir.DEG, suffix.deg)

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
cpms <- cpms.median ## Choose one normalization factor for the CPMs used below
cpms.log10 <- log10(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts

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
pdf(file=paste(sep="", all.prefix,"_sample_correl_heatmap_counts.pdf"))
hm <- heatmap.2(as.matrix(cor(all.counts.mapped)),  scale="none", trace="none", 
                main="Correlation between raw counts",
                col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
quiet <- dev.off()

## Draw sample correlation heatmaps for CPM
#heatmap(cor(all.counts.mapped), scale = "none")
pdf(file=paste(sep="", all.prefix,"_sample_correl_heatmap_cpms.pdf"))
hm <- heatmap.2(as.matrix(cor(cpms)),  scale="none", trace="none", 
                main="Correlation between CPM",
                col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
quiet <- dev.off()

## Plot the first versus second components of samples
cpms.pc <- prcomp(t(cpms))
pdf(file=paste(sep="", all.prefix,"_CPM_PC1-PC2.pdf"))
plot(cpms.pc$x[,1:2], panel.first=grid(), type="n", main="First components from PCA-transformed CPMs")
text(cpms.pc$x[,1:2], labels = sample.conditions, col=cols.samples)
quiet <- dev.off()


################################################################
## Analyse between-replicate reproducibility
################################################################
verbose("Plotting betwen-replicate comparisons", 1)

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
  
  ## Create a specific result directory for this differential analysis
  dir.analysis <- file.path(dir.DEG, paste(sep="", cond1, "_vs_", cond2))
  dir.create(path = dir.analysis, showWarnings = FALSE, recursive = TRUE)
  dir.figures <- file.path(dir.analysis, "figures")
  dir.create(path = dir.figures, showWarnings = FALSE, recursive = TRUE)
  
  ## Only keep genes detected in at least min.rep samples, which is defined as 
  ## the minimal number of replicates per condition.
  min.rep <- min(length(samples1), length(samples2))
  current.counts <- current.counts[rowSums(current.counts > 1) >= min.rep,]
  # dim(current.counts)
  
  ################################################################
  ## Compare counts per million reads (CPMs) between samples. 
  ## This is just to get an intuitive idea, since CPMs are 
  ## not recommended for diffferential detection.
  verbose("\t\tmean CPM plot", 2)
  mean.cpm1 <- apply(cpms[,samples1],1, mean) + epsilon
  mean.cpm2 <- apply(cpms[,samples2],1, mean) + epsilon
  pdf(file=file.path(dir.figures, paste(sep = "", "CPM_plot_", cond1, "_vs_", cond2, ".pdf")))
  plot(mean.cpm1, mean.cpm2, 
       main = "Mean counts per million reads (log scales)",
       xlab=paste(cond1),
       ylab=paste(cond2),
       col="darkblue",log="xy",
       panel.first=grid(lty="solid", col="#DDDDDD"))
  abline(a=0, b=1)
  dev.off()
  
  ## Draw MA plot with CPMs
  verbose("\t\tCPM MA plot", 2)
  A <- log2(mean.cpm1*mean.cpm2)/2
  M <- log2(mean.cpm1/mean.cpm2)
  pdf(file=file.path(dir.figures, paste(sep = "", "CPM_MA_plot_", cond1, "_vs_", cond2, ".pdf")))
  plot(A, M,  
       main = paste(cond1, "vs", cond2, ": CPMs MA plot"),
       xlab=paste(sep="", "A = log2(", cond1, "*", cond2, ")/2"),
       ylab=paste(sep="", "M = log2(", cond1, "/", cond2, ")"),
       col="purple",
       panel.first=grid(lty="solid", col="#DDDDDD"))
  abline(h=0)
  dev.off()
  
  
  ## Define the names of the columns for significant genes according to different criteria
  padj.is.selected.column <- paste(sep="", "padj_", FDR.threshold)
  FDR.is.selected.column <- paste(sep="", "FDR_", FDR.threshold)
  evalue.is.selected.column <- paste(sep="", "Evalue_", Evalue.threshold)
  
  
  ################################################################
  ## DESeq2 analysis
  ################################################################

  verbose("\t\tDESeq2 analysis", 2)
  
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
  
  
  deseq2.dds <- DESeq(deseq2.dds)      ## Differential analysis with negbin distrib
  deseq2.res <- results(deseq2.dds, independentFiltering=FALSE, pAdjustMethod = "BH")  ## Collect the result table
  deseq2.res.sorted <- deseq2.res[order(deseq2.res$padj),]    ## Sort result table by adjusted p-value
  # head(deseq2.res.sorted)
  # tail(deseq2.res.sorted)
  # summary(deseq2.res.sorted)
  
  ## Compute E-value
  deseq2.res.sorted$Evalue <- deseq2.res.sorted$pvalue * nrow(deseq2.res.sorted)
  mcols(deseq2.res.sorted)[7,"type"] <- "results"
  mcols(deseq2.res.sorted)[7,"description"] <- "Expected nb of FP (=pvalue * number of tests)"
  
  ## Compute the rank of p-value sorted genes
  deseq2.res.sorted$pval_rank <- rank(deseq2.res.sorted$pvalue)
  mcols(deseq2.res.sorted)[8,"type"] <- "results"
  mcols(deseq2.res.sorted)[8,"description"] <- "Rank of genes sorted by increasing p-values"
  
  ## Indicate the sense of the regulation (up-down)
  deseq2.res.sorted$sign <- sign(deseq2.res.sorted$log2FoldChange)
  mcols(deseq2.res.sorted)[9,"type"] <- "results"
  mcols(deseq2.res.sorted)[9,"description"] <- "Sign of the regulation (1= up, -1 = down)"
  
  
  ## Label the genes passing the FDR and E-value thresholds
  deseq2.res.sorted[[padj.is.selected.column]] <- (deseq2.res.sorted$padj < FDR.threshold)*1
  mcols(deseq2.res.sorted)[10,"type"] <- "results"
  mcols(deseq2.res.sorted)[10,"description"] <- paste("padj <", FDR.threshold)

  deseq2.res.sorted[[evalue.is.selected.column]] <- (deseq2.res.sorted$Evalue < Evalue.threshold)*1
  mcols(deseq2.res.sorted)[11,"type"] <- "results"
  mcols(deseq2.res.sorted)[11,"description"] <- paste("Evalue <", Evalue.threshold)
  
  mcols(deseq2.res.sorted,use.names=TRUE)
  
  ## Generate and export a result table.
  ## - round numerical values to 3 significant digits.
  ## - add a column with gene_ids, to enable exporting a column header.
  DESeq2.result.table <- cbind("gene_id" = row.names(deseq2.res.sorted),
                               data.frame(deseq2.res.sorted))
  for (col in setdiff(names(DESeq2.result.table), c("gene_id"))) {
    DESeq2.result.table[,col] <- signif(digits=3, DESeq2.result.table[,col])
  }
  deseq2.result.file <- file.path(dir.analysis, 
                        paste(sep = "", cond1, "_vs_", cond2, 
                              "_", suffix.DESeq2, ".tab"))
  write.table(x = DESeq2.result.table, row.names = FALSE,
              file = deseq2.result.file, sep = "\t", quote=FALSE)
  
  verbose(paste(sep="", "\t\tDESeq2 result file\t", deseq2.result.file), 1)
  
  ################################################################
  ## Export DESeq2 plots
  verbose("\t\tExporting DESeq2 plots", 2)
  
  ## Histogram of the nominal p-values
  pdf(file=file.path(dir.figures, paste(sep = "", "DESeq2_pval_hist_", cond1, "_vs_", cond2, ".pdf")))
  hist(deseq2.res.sorted$pvalue, breaks=seq(from=0, to=1, by=epsilon),
       xlab="Nominal p-value",
       ylab="Number of genes",
       main=paste(cond1, "vs", cond2, "; DESeq2 p-values"), col="purple")
  quiet <- dev.off()
  
  
  # Transform the raw discretely distributed counts to apply 
  # some multivariate methods such as PCA, clustering etc.
  deseq2.rld <- rlogTransformation(deseq2.dds, blind=TRUE)

  ## MA plot
  verbose("\t\t\tDESeq2 MA plot", 2)
  MA.ymax <- min(4,max(deseq2.res$log2FoldChange))
  MA.ymin <- max(-4,min(deseq2.res$log2FoldChange))
  MA.ylim <- c(MA.ymin, MA.ymax)
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_plotBCV_", cond1, "_vs_", cond2, "_MAplot.pdf")))
  plotMA(deseq2.dds, main="DESeq2", ylim=MA.ylim)
  grid(lty="solid", col="#CCCCCC")
  quiet <- dev.off()
  
  # PCA plot of the samples
  verbose("\t\t\tDESeq2 PCA plot", 2)
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_plotBCV_", cond1, "_vs_", cond2, "_pca_DESeq2.pdf")))
  print(plotPCA(deseq2.rld, intgroup=c("condition")))
  quiet <- dev.off()
  
  ################################################################
  ## edgeR analysis
  ################################################################
  verbose("\t\tedgeR analysis", 2)
  
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
  
  ## Compute the E-value
  edger.tt$table$Evalue <- edger.tt$table$PValue * nrow(edger.tt$table)
  
  ## Compute the rank of p-value sorted genes
  edger.tt$table$pval_rank <- rank(edger.tt$table$PValue)
  
  ## Indicate the sense of the regulation (up-down)
  edger.tt$table$sign <- sign(edger.tt$table$logFC)
  
  ## Label the genes passing the FDR and E-value thresholds
  edger.tt$table[, FDR.is.selected.column] <- (edger.tt$table$FDR < FDR.threshold)*1
  edger.tt$table[, evalue.is.selected.column] <- (edger.tt$table$Evalue < Evalue.threshold)*1
  
  #  head(edger.tt) 
  
  
  ## Generate and export a result table.
  ## - round numerical values to 3 significant digits.
  ## - add a column with gene_ids, to enable exporting a column header.
  edgeR.result.table <- cbind(data.frame("gene_ID"= row.names(edger.tt$table)), 
                        signif(edger.tt$table, digits=3))
  
  edgeR.result.file <- file.path(dir.analysis, 
                        paste(sep = "", cond1, "_vs_", cond2, 
                              "_", suffix.edgeR, ".tab"))
  write.table(x = edgeR.result.table, row.names = FALSE,
              file = edgeR.result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tedgeR result file\t", edgeR.result.file), 1)
  
  ################################################################
  ## Export edgeR plots
  
  verbose("\t\tGenerating edgeR plots", 2)
  
  ## Histogram of the nominal p-values
  verbose("\t\t\tedgeR nominal p-values histogram", 2)
  pdf(file=file.path(dir.figures, paste(sep = "", "edgeR_pval_hist_", cond1, "_vs_", cond2, ".pdf")))
  hist(edger.tt$table$PValue, breaks=seq(from=0, to=1, by=epsilon),
       xlab="Nominal p-value",
       ylab="Number of genes",
       main=paste(cond1, "vs", cond2, "; edgeR p-values"), col="purple")
  quiet <- dev.off()
  
  ## Smear plot
  verbose("\t\t\tedgeR Smear plot", 3)
  edger.deg <- rownames(edger.tt$table)[edger.tt$table$FDR < FDR.threshold] ## List of differentially expressed genes reported by edgeR
  pdf(file=file.path(dir.figures, paste(sep = "", "edgeR_plotSmear_", cond1, "_vs_", cond2, ".pdf")))
  plotSmear(d, de.tags = edger.deg)
  quiet <- dev.off()
  
  ## MDS plot (Multidimensional scaling plot of distances between gene expression profiles)
  verbose("\t\t\tedgeR MDS plot", 3)
  pdf(file=file.path(dir.figures, paste(sep="", "edgeR_MDS_plot_", cond1, "_vs_", cond2, ".pdf")))
  plotMDS(d, labels=current.labels, 
          col=c("darkgreen","blue")[factor(sample.conditions[names(current.counts)])])
  quiet <- dev.off()
  
  # Mean-variance relationship
  verbose("\t\t\tedgeR mean-variance relationship", 3)
  pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotMeanVar_", cond1, "_vs_", cond2, ".pdf")))
  plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
  quiet <- dev.off()
  
  
  # BCV (Biological Coefficient of Variation)
  verbose("\t\t\tedgeR BCV plot", 3)
  pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotBCV_", cond1, "_vs_", cond2, ".pdf")))
  plotBCV(d)
  quiet <- dev.off()
  
  ################################################################
  ## Merge DESeq2 and edgeR result tables + add some row-wise statistics. 
  ## Beware: gene orders are not the same since they are sorted by adjusted p-value. 
  verbose("\t\tExporting merged table with DESeq2 and edgeR results", 2)
  gene.ids <- row.names(current.counts) ## Use the same order for gene IDs as in the original count table
  result.table <- cbind("gene_id" = DESeq2.result.table[gene.ids,1:1], 
                        "n_min" = apply(current.counts, 1, min),
                        "n_mean" = apply(current.counts, 1, mean),
                        "n_median" = apply(current.counts, 1, median),
                        "n_max" = apply(current.counts, 1, max),
                        "n_sd" = apply(current.counts, 1, sd),
                        "n_sd_mean_ratio" = apply(current.counts, 1, sd)/apply(current.counts, 1, mean),
                        "n_var" = apply(current.counts, 1, var),
                        "n_var_mean_ratio" = apply(current.counts, 1, var)/apply(current.counts, 1, mean),
                        DESeq2 = DESeq2.result.table[gene.ids,2:ncol(DESeq2.result.table)],
                        edgeR = edgeR.result.table[gene.ids,2:ncol(edgeR.result.table)])
  result.table[, paste(sep="", "both_padj_", FDR.threshold)] <- 
    1*(result.table[,paste(sep="", "edgeR.FDR_", FDR.threshold)] & result.table[,paste(sep="", "DESeq2.padj_", FDR.threshold)])
  result.table[, paste(sep="", "edgeR_only_padj_", FDR.threshold)] <- 
    1*(result.table[,paste(sep="", "edgeR.FDR_", FDR.threshold)] & !result.table[,paste(sep="", "DESeq2.padj_", FDR.threshold)])
  result.table[, paste(sep="", "DESeq2_only_padj_", FDR.threshold)] <- 
    1*(!result.table[,paste(sep="", "edgeR.FDR_", FDR.threshold)] & result.table[,paste(sep="", "DESeq2.padj_", FDR.threshold)])
  result.table[, paste(sep="", "none_padj_", FDR.threshold)] <- 
    1*(!result.table[,paste(sep="", "edgeR.FDR_", FDR.threshold)] & !result.table[,paste(sep="", "DESeq2.padj_", FDR.threshold)])
  result.table[, evalue.is.selected.column] <- 
    1*(result.table[,paste(sep="", "edgeR.Evalue_", Evalue.threshold)] & result.table[,paste(sep="", "DESeq2.Evalue_", Evalue.threshold)])
  
  result.file <- file.path(dir.analysis, 
                           paste(sep = "", cond1, "_vs_", cond2, 
                                 "_", suffix.deg, "_DESeq2_and_edgeR.tab"))
  write.table(x = result.table, row.names = FALSE,
              file = result.file, sep = "\t", quote=FALSE)
  verbose(paste(sep="", "\t\tMerged result file\t", result.file), 1)
  
  
  ## Compare DESeq2 and edgeR nominal p-values
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_vs_edgeR_pvalues_", cond1, "_vs_", cond2, ".pdf")))
  plot(result.table$edgeR.PValue, result.table$DESeq2.pvalue,
       xlab="edgeR nominal p-value (log scale)", ylab="DESeq2 nominal p-value (log scale)",
       log="xy", main=paste(cond1, "vs", cond2, "; P-value comparisons"),
       col="#888888",
       panel.first=grid(lty="solid", col="#DDDDDD"))
  abline(a=0, b=1)
  #  abline(h=FDR.threshold)
  #  abline(v=FDR.threshold)
  both <- result.table[, paste(sep="", "both_padj_", FDR.threshold)] == 1
  edgeR.only <- result.table[, paste(sep="", "edgeR_only_padj_", FDR.threshold)] == 1
  DESeq2.only <- result.table[, paste(sep="", "DESeq2_only_padj_", FDR.threshold)] == 1
  points(result.table$edgeR.PValue[both],result.table$DESeq2.pvalue[both], col="darkgreen")
  points(result.table$edgeR.PValue[DESeq2.only],result.table$DESeq2.pvalue[DESeq2.only], col="red")
  points(result.table$edgeR.PValue[edgeR.only],result.table$DESeq2.pvalue[edgeR.only], col="orange")
  legend("topleft", 
         col=c("darkgreen", "red", "orange", "#888888"), pch=1, cex=1, 
         bg="white", bty="o",
         legend = c(paste(sum(both, na.rm=TRUE), "both"),
                    paste(sum(DESeq2.only, na.rm=TRUE), "DESeq2 only"),
                    paste(sum(edgeR.only, na.rm=TRUE), "edgeR only"),
                    paste(sum(result.table[, paste(sep="", "none_padj_", FDR.threshold)], na.rm=TRUE), "none")))
  dev.off()
  
  ## Draw Venn diagram with number of genes declared significant with edgeR and DESeq2, resp
  venn.counts <- vennCounts(result.table[,c(paste(sep="", "edgeR.FDR_", FDR.threshold), 
                                            paste(sep="", "DESeq2.padj_", FDR.threshold))])
  
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_vs_edgeR_Venn_", cond1, "_vs_", cond2, ".pdf")))
  vennDiagram(venn.counts, cex=1, main=paste(cond1, "vs", cond2, "; FDR <", FDR.threshold))  
  dev.off()
  
  ################################################################
  ## Functional enrichment analysis
  gene.IDs.edgeR.FDR <- gene.ids[edgeR.result.table[,FDR.is.selected.column]==1]
#  library("clusterProfiler")
#  library("org.EcK12.eg.db")
#  eg = bitr(gene.IDs.edgeR.FDR, fromType="SYMBOL", toType="ENTREZID", annoDb="org.EcK12.eg.db")
  #eg = bitr(gene.IDs.edgeR.FDR, fromType="SYMBOL", toType="ENTREZID", annoDb="ecoliK12.db0")
  
  ################################################################
  ## Summarise results of the current analysis
  
  ## Instantiate a data frame for the current analysis
  current.summary <- data.frame(
    "analysis"=paste(sep="", cond1, "_vs_", cond2),
    "cond1" = cond1, "cond2" = cond2)

  ## DESeq2 results
  current.summary[,paste(sep="", "DESeq2.padj_", FDR.threshold)] <- 
    sum(DESeq2.result.table[,padj.is.selected.column], na.rm=TRUE)
  current.summary[,paste(sep="", "edgeR.FDR_", FDR.threshold)] <-
    sum(edgeR.result.table[,FDR.is.selected.column])
  current.summary[,paste(sep="", "DESeq2.Evalue_", Evalue.threshold)] <-
    sum(DESeq2.result.table[,evalue.is.selected.column], na.rm=TRUE)
  current.summary[,paste(sep="", "edgeR.Evalue_", Evalue.threshold)] <-
    sum( edgeR.result.table[,evalue.is.selected.column])
  
  if (i == 1) {
    summary.per.analysis <- current.summary
  } else {
    summary.per.analysis <- rbind(summary.per.analysis, current.summary)
  }
  

}

## Export summary table
verbose("Exporting summary table", 1)
summary.file <- paste(sep="", all.prefix, "_summary_per_analysis.tab")
write.table(x = summary.per.analysis, row.names = FALSE,
            file = summary.file, sep = "\t", quote=FALSE)
verbose(paste(sep="", "\tSummary per analysis\t", summary.file), 1)


