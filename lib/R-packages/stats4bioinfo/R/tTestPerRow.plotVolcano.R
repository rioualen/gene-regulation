#' @title Volcano plot from a tTestPerRow result object.
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Draw a volcano plot from the table produced by stats4bioinfo::tTestPerRow().
#' Abcissa represents the mean difference between groups.
#' Ordinate the significance of the Welch test (sig=-log10(e-value)).
#' @param ttpr.result Must be the result obtained with tTestPerRow()
#' @param control.type="fdr"   Type of control (supported: "fdr", "e.value", "p.value")
#' @param alpha=0.05    Alpha threshold for the control of false positives
#' @param ... Additional parameters are passed to plot()
#' @return no return object
#' @examples
#' ## Parameters to generate and analyse the simulated dataset
#' n.h1 <- 100 ## Number of rows under H1
#' n.h0 <- 100 ## Number of rows under H0
#' alpha <- 0.05
#' sample.labels <- c(rep("a", 50), rep("b", 50))
#' 
#' ## Generate an artificial dataset with n.h1 truly positive and n.h0 truly negative rows
#' samples.per.groups <- unlist(table(sample.labels))
#' h1.data <- rnormPerGroup(n=samples.per.groups, mean=c(0, 0.5), sd=c(1,1), nrow=n.h1)
#' h0.data <- rnormPerGroup(n=samples.per.groups, mean=c(0, 0), sd=c(1,1), nrow=n.h0)
#' x <- rbind(h1.data$x, h0.data$x)
#' x.status <- rep(c("H1","H0"), times=c(n.h1, n.h0))
#' table(x.status)
#' 
#' ## Run Student t-test on the simulated data set
#' x.student <- tTestPerRow(x, var.equal=TRUE, cl=sample.labels, alpha=alpha, m2.minus.m1=TRUE)
#' 
#' ## Draw a classical Volcano plot with -log10(p-value) on the Y axis. 
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", control.type = "p.value")
#' 
#' ## Add a threshold on the effect size
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", control.type = "p.value", effect.threshold = 0.6)
#' 
#' ## Draw a Volcano plot with horizontal bars denoting the confidence intervals 
#' ## around the difference between means
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", 
#'    plot.ci=TRUE, control.type="p.value", alpha=)
#' 
#' ## Draw a volcano plot where the colors represent the true status (TP, FP, TN, FN) 
#' ## rather than the significance.
#' (confusion.table <- table(x.status, positive=x.student$table$p.value < alpha))
#' pred.status <- rep(NA, times=nrow(x.student$table))
#' pred.status[x.status=="H1" & x.student$table$p.value < alpha] <- "TP" ## True positives
#' pred.status[x.status=="H0" & x.student$table$p.value >= alpha] <- "TN" ## True negatives
#' pred.status[x.status=="H0" & x.student$table$p.value < alpha] <- "FP" ## False positive
#' pred.status[x.status=="H1" & x.student$table$p.value >= alpha] <- "FN" ## False negatives
#' table(pred.status)
#' pred.status.colors <- c("TP"="darkgreen", "FP"="red", "FN"="orange", "TN"="grey")
#' tTestPerRow.plotVolcano(x.student, legend.corner=NULL, 
#'    plot.ci=TRUE, plot.points=TRUE, tick.size=0.08, control.type="p.value", 
#'    col.points=pred.status.colors[pred.status])
#'  legend("topleft", legend=paste(table(pred.status), names(table(pred.status))),
#'         col=pred.status.colors[names(table(pred.status))], lwd=2)
#' 
#' ## Draw two volcano plots separating positive and negative tests, to show that 
#' ## the p-value cutoff is equivalent to a selection of the confidence interval 
#' ## that do not cross the 0 value. Set the p-value threshold to 1/N, which is 
#' ## equivalent to set the control on e-value <= 1.
#' par(mfrow=c(1,2))
#' corrected.alpha <- 1/nrow(x.student$table)
#' tTestPerRow.plotVolcano(x.student, col.points=NULL, col.positive="black",
#'    legend.corner=NULL, plot.ci=TRUE, tick.size=0, 
#'    control.type="p.value", alpha=corrected.alpha, main="Positive tests")
#' tTestPerRow.plotVolcano(x.student, col.points="grey", col.positive=NA,
#'    legend.corner=NULL, plot.ci=TRUE, tick.size=0, 
#'    control.type="p.value", alpha=corrected.alpha, main="Negative tests")
#' par(mfrow=c(1,1))
#' 
#' @export
tTestPerRow.plotVolcano <- function(
  ttpr.result,
  control.type="fdr",
  alpha=ttpr.result$alpha,
  ylab=paste(sep="", "-log10(", control.type, ")"),
  plot.ci=FALSE,
  ... ## additional parameters are passed to the VolcanoPlot function
) {
  
  multi.t <- ttpr.result$table
  
  ## Identify the positive tests
  if (alpha != ttpr.result$alpha) {
    warning(paste("Selecting positive with alpha=", alpha, "different from ttpr alpha=", ttpr.result$alpha))
    if (plot.ci) {
      warning(paste(sep="", "Computing ",100*(1-alpha),"% confidence intervals"))
      ## Recompute the confidence interval width using user-provided alpha
      ## Compute the width of a confidence interval around the difference between means
      if (ttpr.result$alternative == "two.sided") {
        critical.t <- qt(p=1-alpha/2, df=multi.t$df)
      } else {
        critical.t <- qt(p=1-alpha, df=multi.t$df)
      }
      multi.t[, "ci.width"] <- 2*critical.t * multi.t$diff.stder
    }
  }
  
  VolcanoPlot(multitest.table=multi.t, control.type=control.type, alpha=alpha, effect.size.col="means.diff", ylab=ylab, plot.ci=plot.ci, ...)
}

#' @title Draw a Volcano plot.
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Draw a volcano plot from a table containing at least one column for the effect size, and another one for the p-value or an equivalent measure of significance (FDR, e-value, FWER, ...).
#' @param result.table A data frame containing one row per feature, and one column per statistics.
#' @param effect.size.col A column number or name, inidicating which column of the result table contains the effect size, which will be 
#' @param control.type="p.value"  A column number or name, indicating which column of the result table contains the p-value or an equivalent indication of the significance of each feature (example: "fdr", "e.value", "p.value")
#' @param alpha=0.05    Alpha threshold for the control of false positives
#' @param effect.threshold=NULL Threshold on the absolute value of the effect size.
#' @param sort.by.pval=FALSE Sort row by p-value in order to plot significant elements on top of non-significant
#' @param xlab="Effect size" Label for the X axis.
#' @param ylab="sig = -log10(p-value)" Label for the Y axis
#' @param xlim    Range of the X axis.
#' @param ylim    Range of the Y axis.
#' @param col.points='#888888' Color(s) for the points. can be either a single value (same color for all points), or a vector of the same length as the numer of genes (rows) of the input table. 
#' @param col.positive='#008800' Color to highlight significant points (called positive). When NULL, positive points are not displayed. 
#' @param col.lines='blue'  Color for the line highlighting the significance threshold
#' @param col.grid=#CCCCCC   Grid color
#' @param pch.positive=3   Point shape for the positive tests (default: 3, i.e. + symbol)
#' @param pch.negative=1   Point shape for the negative tests (default: 1, i.e. o symbol)
#' @param legend.corner="bottomleft"  Corner for the legend. When NULL, no legend is displayed.
#' @param full.legend=TRUE Plot additional indications on the legend (number of elements passing the different adjusted p-values)
#' @param legend.cex=1   Font size for the legend.
#' @param cex=0.6     Point size, passed to plot() and lines()
#' @param plot.points=TRUE Plot one point per row with the significance (sig=-log10(p-value)) as a function of the effect size.
#' @param plot.ci=FALSE Plot horizontal bars to denote confidence intervals around the mean difference.
#' @param tick.size=0.05 Height of the vertical lines denoting the boundaries of confidence intervals.
#' @param ... Additional parameters are passed to plot()
#' @return no return object
#' 
#' @export
VolcanoPlot <- function(
  multitest.table,
  effect.size.col,
  control.type="p.value",
  alpha=0.05,
  effect.threshold=NULL,
  sort.by.pval=FALSE,
  plot.points=TRUE,
  plot.ci=FALSE,
  xlab="Effect size",
  ylab=paste(sep="", "-log10(",control.type,")"),
  col.points = '#888888',
  col.positive='#008800',
  col.grid = '#BBBBBB',
  col.lines = 'blue',
  pch.positive=3,
  pch.negative=1,
  xlim = NULL,
  ylim = NULL,
  cex=0.6,
  legend.corner="bottomleft",
  full.legend=FALSE,
  legend.cex=1,
  tick.size=0.05,
  #                         Y.score = "sig", ## Score to plot on the Y axis. Supported: sig (default), p.value e.value
  ... ## additional parameters are passed to the plot function
) {
  
  ## Identify features declared positive according to the specified alpha
  positive <- multitest.table[,control.type] <= alpha
  # table(positive)
  
  ## Apply threshold on effect size if required
  if (!is.null(effect.threshold)) {
    positive[abs(multitest.table[,effect.size.col]) < effect.threshold] <- FALSE
  }
  
  ## Define point colors and shapes. 
  ## Note: the attributes col.points and pch.points can either be either 
  ## a single value (for all points) or a  vector with one user-specified 
  ## color per point.
  multitest.table$color <- col.points
  if (length(col.points) != nrow(multitest.table)) {
    if ((is.na(col.positive)) | (is.null(col.positive))) {
      multitest.table[positive, "color"] <- rep(NA, times=sum(positive))
    } else {
      multitest.table[positive, "color"] <- col.positive
    }
  }
  multitest.table$pch <- pch.negative
  if (!is.null(pch.positive)) {
    multitest.table[positive, "pch"] <- pch.positive
  }
  
  ## Sort the table before plotting, to ensure that positive points appear visible on top of negative points
  if (sort.by.pval) {
    order <- order(multitest.table[,control.type], decreasing=TRUE)
    multitest.table.sorted <- multitest.table[order,]
    if (length(col.points) == nrow(multitest.table)) {
      col.points <- col.points[order]
    }
  } else {
    multitest.table.sorted <- multitest.table
  }
  
  ## Select Y values depending on the control type (p.value, e.value, fdr)
  y.values.ori <- -log10(multitest.table.sorted[, control.type])
  
  ## Fix a problem with infinite Y values resulting from 0 values for the control (p-value or derived stat)
  y.values <- y.values.ori
  y.value.max <- 320 ## Maximal value for Y corresponds to the precision of floating point computation for the p-values (~ 1e-320)
  y.values[is.infinite(y.values)] <- y.value.max
  
  ## Replace NA values by -1, so they appear outside of the plot
  y.values[is.na(y.values)] <- -1
  
  ## Compute confidence interval limits
  if (plot.ci) {
    ci.min <- multitest.table.sorted[, effect.size.col] - multitest.table.sorted[, "ci.width"]/2
    ci.max <- multitest.table.sorted[, effect.size.col] + multitest.table.sorted[, "ci.width"]/2
  }  
  
  ## Define limits of X and Y scales
  if (is.null(xlim)) {
    if (plot.ci) {
      xmax <- max(abs(c(ci.min, ci.max)))
    } else {
      xmax <- max(abs(multitest.table[, effect.size.col]))
    }
    xlim<- c(-xmax, xmax)*1.05 ## Add a 5% margin on X axis
  }
  if (is.null(ylim)) {
    ylim <- c(min(y.values), max(1,y.values))
  }
  
  ## Identify the points above Y limits, to denote them by a different symbol
  above.ymax <- y.values.ori > ylim[2]
  y.values[above.ymax] <- ylim[2]
  multitest.table.sorted$pch[above.ymax] <- 17
  multitest.table.sorted$color[above.ymax] <- "purple"
  
  ################################################################
  ## Draw the Volcano plot
  plot(multitest.table.sorted[, effect.size.col],y.values,
       xlab=xlab,
       ylab=ylab,
       cex=cex,
       xlim=xlim,
       ylim=ylim,
       type="n",
       panel.first=grid(lty='solid', col=col.grid),
       ...)
  
  
  ## Plot the extent + boundaries confidence intervals
  if (plot.ci) {
    arrows(x0 = ci.min, 
           y0 = y.values,
           x1 = ci.max,
           y1 = y.values, 
           col=multitest.table.sorted$color, angle=90, length=tick.size/2, code=3)
  } 
  
  ## Plot the points
  if (plot.points) {
    points(x=multitest.table.sorted[, effect.size.col],
           y=y.values, cex=cex,
           col=multitest.table.sorted$color,
           pch=multitest.table.sorted$pch)
  }
  
  
  ## Vertical line to denote the null position (corresponding to no difference between groups)
  abline(v=0,col="black", lty="dashed", lwd=1) 
  
  ## Add vertical lines to denote the threshold on effect size
  if (!is.null(effect.threshold)) {
    abline(v=c(-effect.threshold, effect.threshold),col=col.lines, lwd=1, lty="solid") ## Vertical line to denote the threshold on effect size
  }
  
  ## Add horizontal lines to denote alpha level
  abline(h=-log10(alpha),col=col.lines, lwd=1, lty="dashed") ## Horizontal line to denote the significance threshold
  
  ## Plot the legend
  if (!is.null(legend.corner)) {
    ## Store legend prameters in a table
    legend.table <- data.frame()
    
    ## Legend for alpha threshold
    if (!is.null(alpha)) {
      legend.table <- rbind(legend.table, 
                            data.frame("name"="alpha", 
                                       "legend"=paste(sep="", control.type," < ", signif(digits=4, alpha)),
                                       "lwd"=2, "col"=col.lines, "pch"=-1, "lty"="dashed"))
    }
    
     
    ## Legend for thresholdon effect size
    if(!is.null(effect.threshold)) {
      legend.table <- rbind(legend.table, 
                            data.frame("name"="effect", 
                                       "legend"=paste(sep="", "abs(", effect.size.col,") >= ", signif(digits=4, effect.threshold)),
                                       "lwd"=2, "col"=col.lines, "pch"=-1, "lty"="solid"))
    }

    ## Legend for the points, depending on their status + display options
    legend.table <- rbind(legend.table, 
                          data.frame("name"="positive", 
                                     "legend"=paste(sep="", sum(positive), " positives"), 
                                     "lwd"=1, "col"= col.positive[1], "pch"=-1, "lty"="blank"))
    legend.table <- rbind(legend.table, 
                          data.frame("name"="negative", 
                                     "legend"=paste(sum(!positive), "negatives"), 
                                     "lwd"=1, "col"= col.points[1], "pch"=-1, "lty"="blank"))

    row.names(legend.table) <- legend.table[, "name"]
    if (plot.points) {
      legend.table["positive", "pch"] <- pch.positive
      legend.table["negative", "pch"] <- pch.negative
    }
  
    
    ## Line type for legend elements
    if (plot.ci) {
      legend.table["positive", "lty"] <- "solid"
      legend.table["negative", "lty"] <- "solid"
    }
    

    
    ## Plot a legend with additional info on the number of positives for different multiple testing corrections
    if (full.legend) {
      legend.pch <- c(legend.table$pch, -1,-1,-1,-1)
      legend(legend.corner, 
             c(paste(sep="","N=",nrow(multitest.table)),
               paste(sep="", sum(positive), " positives (",control.type," <= ", alpha, ")"),
               paste(sum(!positive), "negatives"),
               paste(sep="", "P-value <= ", alpha, ": ", sum(multitest.table[,"p.value"] <= alpha)),
               paste(sep="", "FDR <= ", alpha, ": ", sum(multitest.table[,"fdr"] <= alpha)),
               paste(sep="", "E-value <=", alpha, ": ", sum(multitest.table[,"e.value"] <= alpha)),
               paste(sep="", "E-value <= 1: ", sum(multitest.table[,"e.value"] <= 1))
             ),
             pch=legend.pch, 
             cex=legend.cex,
             lwd=legend.table$lwd, 
             lty=c(legend.table$lty, "blank", "blank", "blank", "blank"),
             col=c(legend.table$col, "white","white","white","white"),
             bty="o", bg="white")
    } else {
      legend(legend.corner, 
             #              c(paste(sep="", control.type," = ", alpha),
             #                paste(sep="", sum(positive), " positives"),
             #                paste(sum(!positive), "negatives")
             #              ),
             legend=as.vector(legend.table$legend),
             pch=as.vector(legend.table$pch), 
             cex=legend.cex,
             lwd=as.vector(legend.table$lwd), 
             lty=as.vector(legend.table$lty),
             col=as.vector(legend.table$col),
             bty="o", bg="white")
      
    }
  }
}


