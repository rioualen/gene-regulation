################################################################
## R code to detect differentially expressed genes using the 
## edgeR BioConductor library. 
##
## Authors: Jeanne Chèneby & Justine Long
## Date: June-July 2015
## 
## Revised by Jacques.van-Helden@univ-amu.fr
## Revision date: 2015-08-15
################################################################

library("edgeR", quietly=TRUE)
#library("limma", quietly=TRUE)
library(gplots, warn.conflicts = FALSE, quietly=TRUE) ## Required for heatmaps.2
library(RColorBrewer, quietly=TRUE)

## Define a color palette for heatmaps
cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))


## The only argument is the file containing all the parameters for the analysis
r.params.path <- commandArgs(trailingOnly = FALSE)[6]
## TEMPORARY FOR DEBUGGING: 
## setwd("~/BeatriceRoche/")
## r.params.path <- "results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_params.R"
source(r.params.path)
setwd(dir.main)

## Read the sample description file, which indicates the 
## condition associated to each sample ID.
sample.desc <- read.delim(sample.description.file, sep="\t", 
                          comment=";", header = TRUE, row.names=1)
sample.ids <- row.names(sample.desc)
sample.conditions <- as.vector(sample.desc[,1])
names(sample.conditions) <- sample.ids

## Read the design file, which indicates the anlayses to be done.
## Each row specifies one differential expression analysis, which 
## consists in comparing two conditions. 
design <- read.delim(design.file, sep="\t", 
                          comment=";", header = TRUE, row.names=NULL)

## Prefix for output files concerning the whole count table (all samples together)
all.prefix <- sub(pattern = ".tab", replacement="", all.counts.table)

## Read the count table
all.counts <- read.delim(all.counts.table, row.names=1, sep="\t")
# names(all.counts)

## Statistics on reads that could not be mapped to genes for different reasons (intergenic, ambiguous, not unique, ...)
all.counts.not.mapped <- all.counts[grep(pattern = "^__", x = row.names(all.counts)), ]
# dim(all.counts.not.mapped)
all.counts.mapped <- all.counts[grep(invert=TRUE, pattern = "^__", x = row.names(all.counts)), ]
# dim(all.counts.mapped)

## Compute the counts per million reads 
## (note: this normalization method is questionable)
cpms <- cpm(all.counts.mapped)    ## Counts per million reads

## Draw sample correlation heatmaps for the raw read counts
pdf(file=paste(sep="", all.prefix,"_sample_correl_heatmap_counts.pdf"))
hm <- heatmap.2(as.matrix(cor(all.counts.mapped)),  scale="none", trace="none", 
                main="Correlation between counts",
                col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
quiet <- dev.off()

## Draw sample correlation heatmaps for CPM
#heatmap(cor(all.counts.mapped), scale = "none")
pdf(file=paste(sep="", all.prefix,"_sample_correl_heatmap_cpms.pdf"))
hm <- heatmap.2(as.matrix(cor(cpms)),  scale="none", trace="none", 
                main="Correlation between counts",
                col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
quiet <- dev.off()

## Plot the first versus second components of samples
pc <- prcomp(t(all.counts.mapped))
pdf(file=paste(sep="", all.prefix,"_PC1-PC2.pdf"))
plot(pc$x[,1:2], panel.first=grid(), type="n")
text(pc$x[,1:2], labels = sample.conditions)
quiet <- dev.off()



## Iterate over analyses
for (i in 1:nrow(design)) {
  
  cond1 <- design[i,1]  ## First condition for the current comparison
  samples1 <- sample.ids[sample.conditions == cond1]
    
  cond2 <- design[i,2]  ## Second condition for the current comparison
  samples2 <- sample.ids[sample.conditions == cond2]
  
  current.samples <- c(samples1, samples2)

  ## Select counts for the samples belonging to the two conditions
  counts <- all.counts.mapped[,current.samples]
  # dim(counts)  ## For test
  
  if (sum(!names(counts) %in% sample.ids) > 0) {
    stop("Count table contains column names without ID in sample description file.")
  }
  
  ## Define conditions and labels for the samples of the current analysis
  current.conditions <- sample.conditions[current.samples]
  current.labels <- paste(current.conditions, names(counts), sep="_")
  
  ## Create a specific result directory for this differential analysis
  dir.analysis <- file.path(dir.DEG, paste(sep="", cond1, "_vs_", cond2))
  dir.create(path = dir.analysis, showWarnings = FALSE, recursive = TRUE)
  dir.figures <- file.path(dir.analysis, "figures")
  dir.create(path = dir.figures, showWarnings = FALSE, recursive = TRUE)
  
  ## Only keep genes detected in at least min.rep samples, which is defined as 
  ## the minimal number of replicates per condition.
  min.rep <- min(length(samples1), length(samples2))
  counts <- counts[rowSums(counts > 1) >= min.rep,]
  # dim(counts)
  
  ################################################################
  ## Convert the count table in a DGEList structure and compute its parameters.
  d <- DGEList(counts=counts, group=sample.conditions[names(counts)])
  d <- calcNormFactors(d)
  d <- estimateCommonDisp(d, verbose=FALSE)
  d <- estimateTagwiseDisp(d, verbose=FALSE)
  
  ################################################################
  ## Detect differentially expressed genes by applying the exact 
  ## test from edgeR package
  edger.de <- exactTest(d)
  
  ## Sort genes by p-value
  edger.tt <- topTags(edger.de, n=nrow(d), sort.by = "PValue")
  
  ## Compute the E-value
  edger.tt$table$Evalue <- edger.tt$table$PValue * nrow(edger.tt$table)

  ## Label the genes passing the FDR and E-value thresholds
  edger.tt$table[, paste(sep="", "FDR_", FDR.threshold)] <- (edger.tt$table$FDR < FDR.threshold)*1
  edger.tt$table[, paste(sep="", "Evalue_", Evalue.threshold)] <- (edger.tt$table$Evalue < Evalue.threshold)*1
  
  #  head(edger.tt) 
  
  #  edger.nc <- cpm(d, normalized.lib.sizes=TRUE)
  edger.deg <- rownames(edger.tt$table)[edger.tt$table$FDR < FDR.threshold]
  
  
  ## Generate and export a result table.
  ## - round numerical values to 3 significant digits.
  ## - add a column with gene_ids, to enable exporting a column header.
  edgeR.result.table <- cbind(data.frame("gene_ID"= row.names(edger.tt$table)), 
                        signif(edger.tt$table, digits=3))
  
  deg.file <- file.path(dir.analysis, 
                        paste(sep = "", cond1, "_vs_", cond2, 
                              "_", suffix.edgeR, ".tab"))
  write.table(x = edgeR.result.table, row.names = FALSE,
              file = deg.file, sep = "\t", quote=FALSE)
  
  ## Summarise results of the current analysis  
  current.summary <- data.frame(
    "analysis"=paste(sep="", cond1, "_vs_", cond2),
    "cond1" = cond1, "cond2" = cond2)
  current.summary[,paste(sep="", "edgeR.FDR_", FDR.threshold)] = sum(edgeR.result.table[,paste(sep="", "FDR_", FDR.threshold)])
  current.summary[,paste(sep="", "edgeR.Evalue_", Evalue.threshold)] = sum(edgeR.result.table[,paste(sep="", "Evalue_", Evalue.threshold)])
  
  if (i == 1) {
    summary.per.analysis <- current.summary
  } else {
    summary.per.analysis <- rbind(summary.per.analysis, current.summary)
  }
  
  ## Plot an histogram of the nominal p-values
  pdf(file=file.path(dir.figures, paste(sep = "", "pval_hist_", cond1, "_vs_", cond2, "_edgeR.pdf")))
  hist(edger.tt$table$PValue, breaks=20)
  quiet <- dev.off()
  
  ################################################################
  ## Export edgeR characteristic plots
  
  ## MA plot
  pdf(file=file.path(dir.figures, paste(sep = "", "edgeR_plotSmear_", cond1, "_vs_", cond2, ".pdf")))
  plotSmear(d, de.tags = edger.deg)
  quiet <- dev.off()
  
  ## MDS plot (Multidimensional scaling plot of distances between gene expression profiles)
  pdf(file=file.path(dir.figures, paste(sep="", "edgeR_MDS_plot_", cond1, "_vs_", cond2, ".pdf")))
  plotMDS(d, labels=current.labels, 
          col=c("darkgreen","blue")[factor(sample.conditions[names(counts)])])
  quiet <- dev.off()
  
  # Mean-variance relationship
  pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotMeanVar_", cond1, "_vs_", cond2, ".pdf")))
  plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
  quiet <- dev.off()
  
  
  # BCV (Biological Coefficient of Variation)
  pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotBCV_", cond1, "_vs_", cond2, ".pdf")))
  plotBCV(d)
  quiet <- dev.off()

}

## Export summary table
summary.file <- file.path(dir.DEG, paste(sep="", suffix.deg, "_summary.tab"))
write.table(x = summary.per.analysis, row.names = FALSE,
            file = summary.file, sep = "\t", quote=FALSE)


