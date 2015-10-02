#' @title Volcano box plot from the result of multiple student tests with resampling (bootstrap, jacknife, sub-sampling).
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Draw a volcano box plot from the object returned by 
#' stats4bioinfo::tTestPerRow.bootstrap(). This representation combines the features 
#' of a volcano plot (significance versus effect size) with the box plots, to indicate 
#' the dispersion (inter-quartille range) and extent of the bootstrapped values.
#' The qbcissa represents the mean difference between groups.
#' The ordinate the significance of the Welch test (sig=-log10(e-value)).
#' Each bootstrapped test is represented by a 2-dimensional box plot, where thick bars 
#' denotes the median, the rectangle delimits the inter-quartile ranges, and the lines
#' indicate the extent of the bootstrapped values. 
#' @param ttpr.result Must be the data frame obtained with tTestPerRow.bootstrap()
#' @param control.type="fdr"   Type of control (supported: "fdr", "e.value", "p.value")
#' @param support.quantile=ttpr.result$support.quantile Support quantile, i.e. the minimal fraction of significant bootstrap iterations required to declare a feature positive. By default, this value is read from the input tTestPerRow() result object. However, it can be over-written in order to estimate the impact of the threshold on support quantile on the result, without having to re-run the analysis.
#' @param plot.sig.boxes=TRUE Draw boxplots for the significance (Y axis)
#' @param plot.effect.boxes=TRUE Draw boxplots for the effect size (X axis)
#' @param plot.rectangles=TRUE  Draw inter-quartile rectangles, i.e. rectangles whose corners are defined by the first and third quartiles, on both axes.
#' @param col.boxplots='#888888' Color(s) for the boxplots. can be either a single value (same color for all boxplots), or a vector of the same length as the numer of multiple tests (rows of the result tables). 
#' @param col.positive='#008800' Color to highlight significant points (called positive).
#' @param col.lines='blue'  Color for the line highlighting the significance threshold
#' @param col.grid=#BBBBBB   Grid color
#' @param lty.positive=3   Line types for the positive tests.
#' @param lty.negative=1   Line types for the negative tests.
#' @param lwd.rectangles= 1 Line width for the IQR rectangles. Set to 0 to avoid drawing rectangle.
#' @param lwd.iqr= 3 Line width for the inter-quartile bars (thick bars from Q1 to Q"). Set to 0 to avoid drawing inter-quartile bars
#' @param lwd.extent=1 Line width for the extent lines Set to 0 to avoid drawing median bars
#' @param lwd.support.quantile=1 Line width for the ticks denoting support quantiles, i.e. the significance corresponding to the specified support quantile. 
#' @param support.ticks=0.03 Width of the ticks denoting the support cutoff. Set to 0 to avoid drawing support cutoff.
#' @param xlab="Effect size"
#' @param ylab="-log10(p-value)"
#' @param xlim    Range of the X axis.
#' @param ylim    Range of the Y axis.
#' @param legend.corner="bottomleft"  Corner for the legend. When NULL, no legend is displayed.
#' @param legend.cex=1   Font size for the legend.
#' @param cex     Point size, passed to plot() and lines()
#' @param ... Additional parameters are passed to plot()
#' @return no return object
#' @examples
#' ## Generate a random set with two samples from distinct 
#' ## populations A and B, characterized by and m rows (features),
#' ## among which m1 are under H1 (mA != mB) and m0 under 
#' ## H0 (mA = mB).
#' m1 <- 500 ## Number of features under H_1
#' m0 <- 500 ## Number of features under H_0
#' m <- m0+m1
#' bootstrap.iterations <- 100
#' 
#' effect.size <- 1 ## Difference between means for the features under H1
#' sample.sizes <- c("A"=30,"B"=30)
#' sample.labels <- rep(names(sample.sizes), times=sample.sizes)
#' table(sample.labels)
#' row.means.g1 <- rep(0, m) ## All rows have a null mean for group 1
#' row.means.g2 <- c(rep(0+effect.size, m1), rep(0, m0)) ## group 2 mean differs depending on the row
#' rand.2grp <- data.frame(cbind(
#'   matrix(nrow=m, ncol=sample.sizes[1], rnorm(m=row.means.g1,sd=1,n=sample.sizes[1]*m)),
#'   matrix(nrow=m, ncol=sample.sizes[2], rnorm(m=row.means.g2,sd=1,n=sample.sizes[2]*m))))
#' 
#' ## Apply Student t-test with classical bootstrap
#' ## (draw samples of same size as original groups, with replacement).
#' ## Set the support quantile to 0.9, in order to select only the features 
#' ## passing the alpha threshold in at least 90% of the iteratios.
#' student.bootstrap <- tTestPerRow.bootstrap(x = rand.2grp, cl=sample.labels, 
#'    iterations=bootstrap.iterations, var.equal=TRUE, support.quantile=0.9)
#'    
#' ## Draw the volcano box plot with default parameters
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap)
#' 
#' ## Draw the volcano box plot with max stringency on support: 
#' ## only consider as positive the rows declared positive in all bootstrap iterations 
#' student.bootstrap.100pcsupport <- tTestPerRow.bootstrap(x = rand.2grp, cl=sample.labels, 
#'    iterations=bootstrap.iterations, var.equal=TRUE, support.quantile=1)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap.100pcsupport)
#' 
#' ## Only draw extent lines and median bars (from quartile 1 to quartile 3) for effect size (horizontal bars and lines)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, lwd.extent=1, lwd.rectangles=0, lwd.iqr=3, 
#'    plot.effect.boxes=TRUE, plot.sig.boxes=FALSE)
#' 
#' ## Only draw extent lines for effect size
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, lwd.extent=1, lwd.rectangles=0, lwd.iqr=0,
#'     plot.effect.boxes=TRUE, plot.sig.boxes=FALSE)
#' 
#' ## Only draw median bars, for both effect size (X axis) and significance (Y axis)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, lwd.extent=0, lwd.rectangles=0, lwd.iqr=3)
#' 
#' ## "Significance bootstrap volano plot": a volcano plot indicating the fluctuations of the significance across the bootstrap iterations
#' ## Compare bootstrap volcanos with different control types
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, lwd.extent=1, lwd.rectangles=0, lwd.iqr=3, support.ticks=0.02,
#'     control.type="p.value", plot.effect.boxes=FALSE, plot.sig.boxes=TRUE, plot.rectangles=FALSE, xlim=c(-0.7,1.7))
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, lwd.extent=1, lwd.rectangles=0, lwd.iqr=3, support.ticks=0.02,
#'     control.type="fdr", plot.effect.boxes=FALSE, plot.sig.boxes=TRUE, plot.rectangles=FALSE, xlim=c(-0.7,1.7))
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, lwd.extent=1, lwd.rectangles=0, lwd.iqr=3, support.ticks=0.02,
#'     control.type="e.value", plot.effect.boxes=FALSE, plot.sig.boxes=TRUE, plot.rectangles=FALSE, xlim=c(-0.7,1.7))
#' 
#' ## Show the options to select subset of graphical elements (rectangles, box plots on each axis). 
#' ## For the sake of comparison, first plot Student t-test on the full dataset (not bootstrapped) + a Volcano CI plot (confidence intervals).
#' ## which is included in the result.
#' x11(width=8, height=12)
#' par(mfrow=c(3,2))
#' par(mar=c(4,4,2,0.5))
#' tTestPerRow.plotVolcano(student.bootstrap$full.set.test, control.type="p.value", legend.corner="topleft", main="Volcano plot")
#' tTestPerRow.plotVolcano(student.bootstrap$full.set.test, control.type="p.value", legend.corner="topleft", plot.ci=TRUE, main="Confidence volcano")
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, control.type="p.value", legend.corner="topleft",
#'    main="Effect size bootstrap volcano", plot.effect.boxes=TRUE, plot.sig.boxes=FALSE, plot.rectangles=FALSE)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, control.type="p.value", legend.corner="topleft",
#'    main="Significance bootstrap volcano", plot.effect.boxes=FALSE, plot.sig.boxes=TRUE, plot.rectangles=FALSE)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, control.type="p.value", legend.corner="topleft",
#'    main="Bootstrap volcano: quartile rectangles", plot.effect.boxes=FALSE, plot.sig.boxes=FALSE, plot.rectangles=TRUE)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, control.type="p.value", legend.corner="topleft",
#'    main="Bootstrap volcano: box plots", plot.effect.boxes=TRUE, plot.sig.boxes=TRUE, plot.rectangles=TRUE)
#' par(mfrow=c(1,1))
#' 
#' 
#' @export
tTestPerRow.bootstrap.VolcanoBoxPlot <- function(
  ttpr.result,
  control.type="fdr",
  support.quantile=ttpr.result$support.quantile,
  plot.sig.boxes=TRUE,
  plot.effect.boxes=TRUE,
  plot.rectangles=TRUE,
  cex=0.7,
  xlab="Effect size",
  ylab=paste(sep="", "-log10(", control.type, ")"),
  # xlab=expression(hat(delta) == bar(X)[1]-bar(X)[2]),
  col.boxplots = '#888888',
  col.positive='#008800',
  col.grid = '#BBBBBB',
  col.lines = 'blue',
  lty.positive="solid",
  lty.negative="solid",
  lwd.rectangles= 1,
  lwd.iqr=3,
  lwd.extent=1,
  lwd.support.quantile=1,
  support.ticks=0.03, 
  xlim = NULL,
  ylim = NULL, 
  legend.corner="topleft",
  legend.cex=0.8,
  #                         Y.score = "sig", ## Score to plot on the Y axis. Supported: sig (default), p.value e.value
  ... ## additional parameters are passed to the plot function
) {
  
  
  ## Check validity of control type
  if (!(control.type %in% c("p.value", "e.value", "fdr"))) {
    stop("Invalid control type for tTestPerRow.bootstrap.VolcanoBoxPlot(). Supported: p.value, e.value, fdr.")
  }
  
  ## Threshold on the number of supporting iterations 
  ## (number of bootstrap tests returning positive result) for the selected control type.
  ## Recompute it because support.quantile can be modified by users for the Volcano plot.
  support.threshold <- ceil(ttpr.result$iterations*support.quantile) 
  
  ## Compute -log10 for all the p-values and derived statistics
  y.columns <- paste(sep="", control.type, c(".min", ".q1", ".median", ".q3", ".max", ".support.quantile"))
  # print(y.columns)
  y.values <- -log10(ttpr.result$stats.per.row[, y.columns])
  names(y.values) <- paste(sep=".", "mlog10", y.columns)
  # View(y.values)
  # head(y.values)
  # dim(y.values)

  ## Compute limits of Y axis
  if (is.null(ylim)) {
    ymax <- max(c(1, unlist(y.values)))
    ymin <- min(unlist(y.values))
    ylim <- c(ymin, ymax)
  }
  
  ## Collect columns for the X axis
  x.columns <- c("effect.size.min","effect.size.q1", "effect.size.median", "effect.size.q3", "effect.size.max")
  x.values <- ttpr.result$stats.per.row[, x.columns]
  # head(x.values)
  # dim(x.values)
  
  ## X limits
  if (is.null(xlim)) {
    xmax <- max(abs(c(x.values$effect.size.min, x.values$effect.size.max)))
    xlim <- c(-xmax, xmax)
  }
  
  ## Identify the positive tests
  control.column <- paste(sep="", control.type, ".support")
  positive <- ttpr.result$stats.per.row[,control.column] >= support.threshold
  # table(positive)
  order <- order(ttpr.result$stats.per.row[,control.column], decreasing=TRUE)

  ## Define point colors and shapes. 
  ## Note: the attributes col.boxplots and lty.points can either be either 
  ## a single value (for all points) or a  vector with one user-specified 
  ## color per point.
  ttpr.result$stats.per.row$color <- col.boxplots
  if (length(col.boxplots) != nrow(ttpr.result$stats.per.row)) {
    if (!is.null(col.positive)) {
      ttpr.result$stats.per.row[positive, "color"] <- col.positive
    }
  }
  # table(ttpr.result$stats.per.row$color)
  ttpr.result$stats.per.row$lty <- lty.negative
  if (!is.null(lty.positive)) {
    ttpr.result$stats.per.row[positive, "lty"] <- lty.positive
  }
  # table(ttpr.result$stats.per.row$lty)
  
  ## Sort the table befor eplotting, to ensure that 
  ## positive points appear visible on top of negative points
  stats.per.row.sorted <- ttpr.result$stats.per.row[order,]
  # stats.per.row.sorted <- ttpr.result$stats.per.row
  
  ################################################################
  ## Draw the Volcano plot
  plot(unlist(x.values$effect.size.median),
       unlist(y.values[, paste(sep="", "mlog10.",control.type, ".median")]),
       xlab=xlab,
       ylab=ylab,
       cex=cex,
       xlim=xlim,
       ylim=ylim,
       type="n",
       panel.first=grid(lty='solid', col=col.grid),
       ...)
  
  
  ## Plot quartile rectangles
  if (plot.rectangles) {
    rect(xleft=x.values$effect.size.q1, 
         ybottom=y.values[, paste(sep="", "mlog10.",control.type, ".q3")], 
         xright=x.values$effect.size.q3, 
         ytop=y.values[, paste(sep="", "mlog10.",control.type, ".q1")], 
         col = NULL, 
         lwd=lwd.rectangles,
         border = ttpr.result$stats.per.row$color, 
         lty = ttpr.result$stats.per.row$lty)
  }  
  
  if (plot.sig.boxes) {
    ## Plot extent lines for p-values
    arrows(x0 = x.values$effect.size.median, 
           y0 = y.values[, paste(sep="", "mlog10.",control.type, ".min")],
           x1 = x.values$effect.size.median,
           y1 = y.values[, paste(sep="", "mlog10.",control.type, ".max")], 
           lwd=lwd.extent,
           lty = ttpr.result$stats.per.row$lty,
           col=ttpr.result$stats.per.row$color, angle=0, length=0, code=3)
    
    ## Plot p-value IQR
    arrows(x0 = x.values$effect.size.median, 
           y0 = y.values[, paste(sep="", "mlog10.",control.type, ".q3")],
           x1 = x.values$effect.size.median,
           y1 = y.values[, paste(sep="", "mlog10.",control.type, ".q1")], 
           lwd=lwd.iqr,
           lty = ttpr.result$stats.per.row$lty,
           col=ttpr.result$stats.per.row$color, angle=0, length=0, code=3)    

    ## Compute positions of support quantile ticks
    if (support.quantile == ttpr.result$support.quantile) {
      quantile.support.ypos <- y.values[, paste(sep="", "mlog10.",control.type, ".support.quantile")]
    } else {
      ## If the requested support quantile differs from the one used during bootstrap, we need to recompute the positions of the support ticks. 
      quantile.support.values <- apply(ttpr.result[[control.type]],1, quantile, support.quantile)
      quantile.support.ypos <- -log10(quantile.support.values)
    }
    ## Plot support quantiles
    arrows(x0 = x.values$effect.size.median - support.ticks/2,
           y0 = quantile.support.ypos,
           x1 = x.values$effect.size.median + support.ticks/2,
           y1 = quantile.support.ypos, 
           lwd=lwd.support.quantile,
           lty = ttpr.result$stats.per.row$lty,
           col=ttpr.result$stats.per.row$color, angle=0, length=0, code=3)    
  }
  
  
  if (plot.effect.boxes) {    
    ## Plot extent lines for effect sizes
    arrows(x0 = x.values$effect.size.min, 
           y0 = y.values[, paste(sep="", "mlog10.",control.type, ".median")],
           x1 = x.values$effect.size.max,
           y1 = y.values[, paste(sep="", "mlog10.",control.type, ".median")], 
           lwd=lwd.extent,
           lty = ttpr.result$stats.per.row$lty,
           col=ttpr.result$stats.per.row$color, angle=0, length=0, code=3)
    
    ## Plot IQR of effect sizes
    arrows(x0 = x.values$effect.size.q1, 
           y0 = y.values[, paste(sep="", "mlog10.",control.type, ".median")],
           x1 = x.values$effect.size.q3,
           y1 = y.values[, paste(sep="", "mlog10.",control.type, ".median")], 
           lwd=lwd.iqr,
           lty = ttpr.result$stats.per.row$lty,
           col=ttpr.result$stats.per.row$color, angle=0, length=0, code=3)
  }    
  
  
  abline(v=0,col="black", lty="dashed")
#  abline(h=0,col="black", lty="dashed")
  abline(h=-log10(ttpr.result$alpha), col=col.lines, lty="solid")
  
  ## Plot the legend
  if (!is.null(legend.corner)) {
    legend.lwd=2
    legend.lty <- c("blank", "blank", lty.positive, lty.negative)
    legend(legend.corner, 
           legend=c(
             paste(sep="", control.type," < ", ttpr.result$alpha),
             paste(sep="", "support >= ", support.threshold, "/", ttpr.result$iterations),
             #             paste(sep="","N=",nrow(ttpr.result$stats.per.row)),
             paste(sep="", sum(positive), " positives"),
             paste(sum(!positive), "negatives")
             ),
           lty=legend.lty, 
           lwd=legend.lwd, 
           col=c("white", "white", col.positive[1], col.boxplots[1]),
           bty="o", bg="white")
  }
}


