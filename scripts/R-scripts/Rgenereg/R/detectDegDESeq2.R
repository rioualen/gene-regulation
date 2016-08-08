################################################################
## DESeq2 analysis
##
## Detect differentially expressed genes (DEG) using the package DESeq2,
## and add a few custom columns (e-value, ...) to the result table.
detectDegDESeq2 <- function(dir.figures=NULL) {
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
  return(deseq2.result.table)
}

