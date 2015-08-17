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

## Define a color palette for heatmaps. I like this Red-Blue palette because 
## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
## - it can be seen by people suffering from red–green color blindness.
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

## Experimental conditions
sample.conditions <- as.vector(sample.desc[,1]) ## Condition associated to each sample
names(sample.conditions) <- sample.ids

## Define a specific color for each distinct condition
exp.conditions <- unique(sample.conditions) ## Set of distinct conditions
cols.conditions <- brewer.pal(length(exp.conditions),"Accent")
names(cols.conditions) <- exp.conditions

## Read the design file, which indicates the anlayses to be done.
## Each row specifies one differential expression analysis, which 
## consists in comparing two conditions. 
design <- read.delim(design.file, sep="\t", 
                          comment=";", header = TRUE, row.names=NULL)

## Prefix for output files concerning the whole count table (all samples together)
## all.prefix <- sub(pattern = ".tab", replacement="", all.counts.table)
all.prefix <- file.path(dir.DEG, suffix.deg)

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


################################################################
## Analyse between-replicate reproducibility
################################################################

for (cond in exp.conditions) {
  ## Create a specific result directory for this condition
  dir.condition <- file.path(dir.DEG, "per_condition", cond)
  dir.create(path = dir.condition, showWarnings = FALSE, recursive = TRUE)
  
  current.samples <- names(all.counts.mapped)[sample.conditions == cond]
  nrep <- length(current.samples)
  max.rep.to.plot <- min(3, nrep)
  
  ## A trick: by adding 0.1 we can view the 0 values on the log plot, at the 0.1 coordinate
  to.plot <- all.counts.mapped[,current.samples][,1:max.rep.to.plot] + 0.1
  
  pdf(file= file.path(dir.condition, paste(sep = "", "between-replicate_compa_plot_", cond1, ".pdf")), width=10, height=10)
  plot(to.plot, log="xy", col=cols.conditions[cond], 
       main=paste(cond, " ; raw counts per replicate (log scale)"))
  dev.off()
}

################################################################
## Run differential expression analysis
################################################################

## Iterate over analyses
for (i in 1:nrow(design)) {
  
  cond1 <- as.vector(design[i,1])  ## First condition for the current comparison
  samples1 <- sample.ids[sample.conditions == cond1]
    
  cond2 <- as.vector(design[i,2])  ## Second condition for the current comparison
  samples2 <- sample.ids[sample.conditions == cond2]
  
  current.samples <- c(samples1, samples2)
  current.counts <- all.counts.mapped[,current.samples]

  ## Select counts for the samples belonging to the two conditions
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
  ## DESeq2 analysis
  ################################################################
  
  ## Create a DESeqDataSet object from the count table + conditions
  condition <- as.factor(as.vector(current.conditions))
  deseq2.dds <- DESeqDataSetFromMatrix(
    countData = current.counts, 
    colData = DataFrame(condition),
    ~ condition)
  
  
  ## Indicate that second condition is the reference condition. 
  ## If not done, the conditions are considered by alphabetical order, 
  ## which may be misleading to interpret the log2 fold changes. 
  deseq2.dds$condition <- relevel(deseq2.dds$condition, cond2) 
  
  
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
  
  ## Label the genes passing the FDR and E-value thresholds
  deseq2.res.sorted[[paste(sep="", "padj_", FDR.threshold)]] <- (deseq2.res.sorted$padj < FDR.threshold)*1
  mcols(deseq2.res.sorted)[8,"type"] <- "results"
  mcols(deseq2.res.sorted)[8,"description"] <- paste("padj <", FDR.threshold)
  deseq2.res.sorted[[paste(sep="", "Evalue_", Evalue.threshold)]] <- (deseq2.res.sorted$Evalue < Evalue.threshold)*1
  mcols(deseq2.res.sorted)[9,"type"] <- "results"
  mcols(deseq2.res.sorted)[9,"description"] <- paste("Evalue <", Evalue.threshold)
  
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
  
  
  ################################################################
  ## Export DESeq2 plots
  
  ## Histogram of the nominal p-values
  pdf(file=file.path(dir.figures, paste(sep = "", "DESeq2_pval_hist_", cond1, "_vs_", cond2, ".pdf")))
  hist(deseq2.res.sorted$pvalue, breaks=seq(from=0, to=1, by=0.01),
       xlab="Nominal p-value",
       ylab="Number of genes",
       main=paste(cond1, "vs", cond2, "; DESeq2 p-values"), col="purple")
  quiet <- dev.off()
  
  
  # Transform the raw discretely distributed counts to apply 
  # some multivariate methods such as PCA, clustering etc.
  deseq2.rld <- rlogTransformation(deseq2.dds, blind=TRUE)

  ## MA plot
  MA.ymax <- min(4,max(deseq2.res$log2FoldChange))
  MA.ymin <- max(-4,min(deseq2.res$log2FoldChange))
  MA.ylim <- c(MA.ymin, MA.ymax)
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_plotBCV_", cond1, "_vs_", cond2, "_MAplot.pdf")))
  plotMA(deseq2.dds, main="DESeq2", ylim=MA.ylim)
  grid(lty="solid", col="#CCCCCC")
  quiet <- dev.off()
  
  # PCA plot of the samples
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_plotBCV_", cond1, "_vs_", cond2, "_pca_DESeq2.pdf")))
  print(plotPCA(deseq2.rld, intgroup=c("condition")))
  quiet <- dev.off()
  
  ################################################################
  ## edgeR analysis
  ################################################################

  ## Convert the count table in a DGEList structure and compute its parameters.
  d <- DGEList(counts=current.counts, group=sample.conditions[names(current.counts)])
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
    
  ################################################################
  ## Export edgeR plots
  
  ## Histogram of the nominal p-values
  pdf(file=file.path(dir.figures, paste(sep = "", "edgeR_pval_hist_", cond1, "_vs_", cond2, ".pdf")))
  hist(edger.tt$table$PValue, breaks=seq(from=0, to=1, by=0.01),
       xlab="Nominal p-value",
       ylab="Number of genes",
       main=paste(cond1, "vs", cond2, "; edgeR p-values"), col="purple")
  quiet <- dev.off()
  
  ## MA plot
  edger.deg <- rownames(edger.tt$table)[edger.tt$table$FDR < FDR.threshold] ## List of differentially expressed genes reported by edgeR
  pdf(file=file.path(dir.figures, paste(sep = "", "edgeR_plotSmear_", cond1, "_vs_", cond2, ".pdf")))
  plotSmear(d, de.tags = edger.deg)
  quiet <- dev.off()
  
  ## MDS plot (Multidimensional scaling plot of distances between gene expression profiles)
  pdf(file=file.path(dir.figures, paste(sep="", "edgeR_MDS_plot_", cond1, "_vs_", cond2, ".pdf")))
  plotMDS(d, labels=current.labels, 
          col=c("darkgreen","blue")[factor(sample.conditions[names(current.counts)])])
  quiet <- dev.off()
  
  # Mean-variance relationship
  pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotMeanVar_", cond1, "_vs_", cond2, ".pdf")))
  plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
  quiet <- dev.off()
  
  
  # BCV (Biological Coefficient of Variation)
  pdf(file= file.path(dir.figures, paste(sep = "", "edgeR_plotBCV_", cond1, "_vs_", cond2, ".pdf")))
  plotBCV(d)
  quiet <- dev.off()
  
  ################################################################
  ## Merge DESeq2 and edgeR result tables + add some row-wise statistics. 
  ## Beware: gene orders are not the same since they are sorted by adjusted p-value. 
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
  result.table[, paste(sep="", "Evalue_", Evalue.threshold)] <- 
    1*(result.table[,paste(sep="", "edgeR.Evalue_", Evalue.threshold)] & result.table[,paste(sep="", "DESeq2.Evalue_", Evalue.threshold)])
  
  result.file <- file.path(dir.analysis, 
                           paste(sep = "", cond1, "_vs_", cond2, 
                                 "_", suffix.deg, "_DESeq2_and_edgeR.tab"))
  write.table(x = result.table, row.names = FALSE,
              file = result.file, sep = "\t", quote=FALSE)
  
  
  ## Compare DESeq2 and edgeR nominal p-values
  pdf(file= file.path(dir.figures, paste(sep = "", "DESeq2_vs_edgeR_pvalues_", cond1, "_vs_", cond2, ".pdf")))
  plot(result.table$edgeR.PValue, result.table$DESeq2.pvalue,
       xlab="edgeR nominal p-value (log scale)", ylab="DESeq2 nominal p-value (log scale)",
       log="xy", main=paste(cond1, "vs", cond2, "; P-value comparisons"),
       col="#888888",
       panel.first=grid(lty="solid", col="#DDDDDD"))
  abline(a=0, b=1)
  abline(h=FDR.threshold)
  abline(v=FDR.threshold)
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
  ## Summarise results of the current analysis
  
  ## Instantiate a data frame for the current analysis
  current.summary <- data.frame(
    "analysis"=paste(sep="", cond1, "_vs_", cond2),
    "cond1" = cond1, "cond2" = cond2)

  ## DESeq2 results
  current.summary[,paste(sep="", "DESeq2.padj_", FDR.threshold)] <- 
    sum(DESeq2.result.table[,paste(sep="", "padj_", FDR.threshold)], na.rm=TRUE)
  current.summary[,paste(sep="", "edgeR.FDR_", FDR.threshold)] <-
    sum(edgeR.result.table[,paste(sep="", "FDR_", FDR.threshold)])
  current.summary[,paste(sep="", "DESeq2.Evalue_", Evalue.threshold)] <-
    sum(DESeq2.result.table[,paste(sep="", "Evalue_", Evalue.threshold)], na.rm=TRUE)
  current.summary[,paste(sep="", "edgeR.Evalue_", Evalue.threshold)] <-
    sum( edgeR.result.table[,paste(sep="", "Evalue_", Evalue.threshold)])
  
  if (i == 1) {
    summary.per.analysis <- current.summary
  } else {
    summary.per.analysis <- rbind(summary.per.analysis, current.summary)
  }
  

}

## Export summary table
summary.file <- paste(sep="", all.prefix, "_summary.tab")
write.table(x = summary.per.analysis, row.names = FALSE,
            file = summary.file, sep = "\t", quote=FALSE)


