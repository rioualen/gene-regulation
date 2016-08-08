detectDEGedgeR <- function(dir.figures=NULL) {
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
    thresholds=thresholds,
    round.digits = 3,
    dir.figures=dir.figures)
  # dim(edger.result.table)
  # dim(edger.result.table)
  # names(edger.result.table)
  # View(edger.result.table)
  return (edger.result.table)
}
