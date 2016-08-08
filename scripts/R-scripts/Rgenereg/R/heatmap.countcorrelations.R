################################################################
## Draw a heatmap with the inter-sample correlation matrix.
heatmap.countCorrelations <- function(count.table,
                                      main="Correlation between raw counts",
                                      plot.file=NULL,
                                      log.transform=FALSE, # Perform a log transformation of the values before plotting
                                      epsilon=0.01, # Add an epsilon to zero values before log transformation, in order to -Inf values
                                      ...
) {

  ## Define a color palette for heatmaps. I like this Red-Blue palette because
  ## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
  ## - it can be seen by people suffering from red–green color blindness.
  cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))


  ## Adapt boxplot size to the number of samples and label sizes
  margin <- max(nchar(names(count.table)))/3+5

  if (log.transform) {
    range(count.table)
    count.table[count.table==0] <- epsilon
    count.table <- log10(count.table)
  }
  count.cor <- as.matrix(cor(count.table))

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }

  hm <- heatmap.2(count.cor,  scale="none", trace="none", breaks=seq(-1,1,length.out = 101),
                  main=main, margins=c(margin,margin),
                  col=cols.heatmap,
                  cellnote = signif(digits=2, count.cor),
                  ...
  )


  if (!is.null(plot.file)) {
    silence <- dev.off()
  }

  return(count.cor)
}

