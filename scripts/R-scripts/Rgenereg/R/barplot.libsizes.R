################################################################
## Draw a barplot with the number of reads per sample
barplot.libsizes <- function(stats.per.sample,
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
  bplt <- barplot(stats.per.sample$Mreads, names.arg = stats.per.sample$label, horiz = TRUE, las=1,
                  xlab="libsum (Million reads per sample)",
                  main="Read library sizes (libsum per sample)",
                  col=stats.per.sample$color)
  grid(col="white", lty="solid",ny = 0)
  text(x=pmax(stats.per.sample$Mreads, 3), labels=stats.per.sample$Mreads, y=bplt,pos=2, font=2)
  if (!is.null(plot.file)) {
    silence <- dev.off()
  }
}

