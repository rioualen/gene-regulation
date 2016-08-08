## Generate a set of plots displaying some sample-wise statistics
sampleDescriptionPlots <- function (sample.desc,
                                    stats.per.sample,
                                    dir.figures,
                                    exploratory.plots=FALSE) {

  plot.files <- c() ## To retun: list of files with plots

  par.ori <- par() ## Save original plot parameters


  ## Library size barplots
  plot.files["Mreads_barplot"] <- file.path(dir.figures, "sample_libsum_barplot.pdf")
  libsize.barplot(stats.per.sample, plot.files["Mreads_barplot"])

  ################################################################
  ## Boxplots of raw counts and CPMs, in linear + log axes.
  ## These plots give a pretty good intuition of the raw data per sample:
  ## library sizes, outliers, dispersion of gene counts.

  ## Boxplot of raw counts
  plot.files["sample_boxplot_counts"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_counts.pdf"))
  count.boxplot(all.counts, stats.per.sample,
                xlab="Raw counts", main="Box plots per sample: raw counts",
                plot.file=plot.files["sample_boxplot_counts"])

  ## Boxplot of log10-transformed counts
  plot.files["sample_boxplot_counts_log10"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_counts_log10.pdf"))
  count.boxplot(all.counts.log10, stats.per.sample,
                xlab="log10(counts)", main="Box plots per sample: log10(counts)",
                plot.file=plot.files["sample_boxplot_counts_log10"])

  ## Boxplot of CPMs
  plot.files["sample_boxplot_CPM"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_CPM.pdf"))
  count.boxplot(cpms, stats.per.sample,
                xlab="CPM", main="Box plots per sample: counts per million reads (CPM)",
                plot.file=plot.files["sample_boxplot_CPM"])

  ## Boxplot of log10-transformed CPMs
  plot.files["sample_boxplot_CPM_log10"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_CPM_log10.pdf"))
  count.boxplot(cpms.log10, stats.per.sample,
                xlab="log10(CPM)", main="Box plots per sample: counts per million reads (CPM)",
                plot.file=plot.files["sample_boxplot_CPM"])

  par <- par.ori ## Restore original plot parameters
  par(mar=c(4.1,5.1,4.1,1.1))

  ## Draw sample correlation heatmaps for the raw read counts
  plot.files["sample_correl_heatmap_counts"] <- paste(sep="", prefix["general.file"],"_sample_correl_heatmap_counts.pdf")
  count.correl.heatmap(all.counts, plot.file=plot.files["sample_correl_heatmap_counts"])
  #   hm <- heatmap.2(,  scale="none", trace="none",
  #                   main="Correlation between raw counts", margins=c(8,8),
  #                   col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))

  ## Draw sample correlation heatmaps for CPM. Actually it gives exactly the
  ## same result as correlation between raw counts, since the correlation has a
  ## standardizing effect.
  # pdf(file=paste(sep="", prefix["general.file"],"_sample_correl_heatmap_cpms.pdf"))
  # hm <- heatmap.2(as.matrix(cor(cpms)),  scale="none", trace="none",
  #                 main="Correlation between CPM",
  #                 col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
  # quiet <- dev.off()

  ## Plot the first versus second components of samples
  cpms.pc <- prcomp(t(cpms))
  plot.file <- paste(sep="", prefix["general.file"],"_CPM_PC1-PC2.pdf")
  plot.files["CPM_PC1-PC2"] <- plot.file
  message("Generating plot", plot.file)
  pdf(file=plot.file)
  plot(cpms.pc$x[,1:2], panel.first=grid(), type="n", main="First components from PCA-transformed CPMs")
  text(cpms.pc$x[,1:2], labels = sample.conditions, col=sample.desc$color)
  quiet <- dev.off()



  ## Exploratory plots, should not be done for all projects.
  if (exploratory.plots) {
    verbose("Drawing generic plots from the whole count table", 1)

    ## Plot the impact of the normalization factor (library sum , median or percentile 75)
    plot.file <- file.path(dir.DEG, paste(sep = "", "CPM_libsum_vs_median_vs_perc75.png"))
    plot.files["CPM_libsum_vs_median_vs_perc75"] <- plot.file
    message("Generating plot", plot.file)
    png(file= plot.file, width=1000, height=1000)
    cols.counts <- as.data.frame(matrix(sample.desc$color, nrow=nrow(all.counts), ncol=ncol(all.counts), byrow = TRUE))
    colnames(cols.counts) <- names(all.counts)
    rownames(cols.counts) <- rownames(all.counts)
    plot(data.frame("libsum" = as.vector(as.matrix(cpms.libsum)),
                    "median" = as.vector(as.matrix(cpms.median)),
                    "perc75" = as.vector(as.matrix(cpms.perc75))),
         col=as.vector(as.matrix(cols.counts)))
    quiet <- dev.off()

    ## Plot some sample-wise statistics
    plot.file <- file.path(dir.DEG, paste(sep = "", "sample_statistics_plots.pdf"))
    plot.files["sample_statistics_plots"] <- plot.file
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=10, height=10)
    par(mar=c(5,5,1,1)) ## adpt axes
    par(mfrow=c(2,2))
    ## Median versus mean
    plot(stats.per.sample[,c("mean", "median")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)

    ## First versus third quartile
    plot(stats.per.sample[,c("perc25", "perc75")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)

    ## Sum versus third quartile.
    plot(stats.per.sample[,c("sum", "perc75")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)

    ## Mean versus third quartile.
    plot(stats.per.sample[,c("mean", "perc75")],
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)
    par(mfrow=c(1,1))
    quiet <- dev.off()
  }

  return(plot.files)
}
