################################################################
## Draw boxplots with read counts per genes for each sample
boxplot.countTable <- function(count.table,
                          sample.desc,
                          xlab="Raw counts",
                          main="Box plots per sample: raw counts",
                          plot.file=NULL) {

  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- max(nchar(sample.desc$label))/3+5
  boxplot.height <- length(sample.ids)/3+2

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }

  par(mar=c(5,boxplot.lmargin,4,1)) ## adapt axes
  boxplot(count.table, horizontal=TRUE, col=sample.desc$color,
          xlab=xlab, names=sample.desc$label,
          main=main, las=1)
  grid(col="grey", lty="solid",ny = 0)
  if (!is.null(plot.file)) {
    silence <- dev.off()
  }
}

