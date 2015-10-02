#' @title Draw an XY plot comparing two series of p-values (or corrected p-values, or e-values).
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Draw an XY plot comparing two series of 
#' p-values (or corrected p-values, or e-values).
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param frame.to.plot       A data.frame comprizing two columns, with the respective p-values (corrected or not) to plot.
#' @param xlab=names(frame.to.plot)[1]  Passed to plot()
#' @param ylab=names(frame.to.plot)[2]  Passed to plot()
#' @param alpha=0.05   Alpha threshold to select positives and compute the contingency table.
#' @param score.name="p-value"  Score name to display on the legend (e.g. "p-value", "e-value", "fdr").
#' @param plot.result=TRUE    Draw a comparison plot
#' @param plot.colors=TRUE    Use colors for the plot   
#' @param plot.cex=0.7     Point size for the plot (passed to graphics::plot() function).
#' @param legend.cex=0.9    Font size factor for legend. This allows to choose separate font sizes to draw the plot itself and the legend.
#' @param legend.corner="bottomright"   Position where legend should be placed.
#' @param ...             Additional parameters are passed to plot()
#'
#' @examples
#' ################################################################
#' ## Load dataset from DenBoer, 2009 and define two groups of interest
#' library(denboer2009)
#' group1 <- "Bo"
#' group2 <- "Bt"
#' 
#' verbose(paste("Selecting samples for groups", group1, "and", group2), 1)
#' selected.groups <- c(group1, group2)
#' selected.samples <- denboer2009.pheno$sample.labels %in% selected.groups
#' 
#' ## Selected expression profiles
#' selected.probesets <- apply(denboer2009.amp == "P", 1, sum) >= 30
#' expr <- denboer2009.expr[selected.probesets, selected.samples]
#' sample.labels <- denboer2009.pheno$sample.labels[selected.samples]
#' verbose(paste("Expression profiles: ", nrow(expr), "probesets x", ncol(expr), "samples"), 1)
#' table(sample.labels)
#' 
#' ################################################################
#' ## Compare p-values returned by Student and Welch t-tests
#' student.result <- tTestPerRow(x=expr,cl = sample.labels, var.equal = TRUE)
#' welch.result <- tTestPerRow(x=expr,cl = sample.labels, var.equal = FALSE)
#' plotPvalCompa(data.frame("Student.e.value"=student.result$table$e.value,
#'                          "Welch.e.value"=welch.result$table$e.value),
#'               score.name = "e-value",
#'               legend.corner="bottomright")
#' 
#' @export
plotPvalCompa <- function (frame.to.plot,
                           xlab=names(frame.to.plot)[1],
                           ylab=names(frame.to.plot)[2],
                           alpha=0.05,
                           score.name="p-value",
                           plot.result=TRUE,
                           plot.colors=TRUE,
                           plot.cex=0.7,
                           legend.cex = 0.9,
                           legend.corner="bottomright",
                           ...) {
  
  
  result <- list()
  
  score1 <- frame.to.plot[,1]
  score2 <- frame.to.plot[,2]
  
  result$score.table <- frame.to.plot
  
  ## Count the number of elements below the threshold in each col
  result$score.table$score1.signif <- score1 <= alpha
  result$score.table$score2.signif <- score2 <= alpha
  result$score.table$score2.or.score1.signif <- result$score.table$score1.signif | result$score.table$score2.signif
  result$score.table$score2.and.score1.signif <- result$score.table$score1.signif & result$score.table$score2.signif
  
  ## Count number of significant rows for each test
  result$nb.signif.score1 <- sum(result$score.table$score1.signif)
  result$nb.signif.score2 <- sum(result$score.table$score2.signif)
  result$nb.signif.both <- sum((result$score.table$score1.signif) & (result$score.table$score2.signif))
  result$nb.signif.any <- sum((result$score.table$score1.signif) | (result$score.table$score2.signif))  
  result$nb.signif.none <- sum((!result$score.table$score1.signif) & (!result$score.table$score2.signif))
  result$jaccard <- result$nb.signif.both / result$nb.signif.any
  
  result$contingency <- table(result$score.table$score1.signif, result$score.table$score2.signif)
  
  ################################################################
  ## Comparison plot
  if (plot.result) {
    
    
    ## Define density colors to highlight places where many points overlap
    density.cols <- densCols(-log10(frame.to.plot), nbin=256,
                             col=colorRampPalette(c("#BBBBBB", "#000000")))
    
    ## Add colors to highlight consistencies and inconsistencies between significant features
    if (plot.colors) {
      density.cols[result$score.table$score1.signif] <- "blue"
      density.cols[result$score.table$score2.signif] <- "#00BB00"
      density.cols[result$score.table$score1.signif & result$score.table$score2.signif] <- "red"
    } else {
      density.cols[result$score.table$score1.signif | result$score.table$score2.signif] <- "#888888"
      density.cols[result$score.table$score1.signif & result$score.table$score2.signif] <- "#444444"
    }
    
    ## Define point characters to highlight consistencies and inconsistencies
    compa.pch <- rep(1, nrow(frame.to.plot))
    compa.pch[result$score.table$score1.signif] <- 3
    compa.pch[result$score.table$score2.signif] <- 4
    compa.pch[result$score.table$score1.signif & result$score.table$score2.signif] <- 8
    
    ## Plot the comparison
    plot(frame.to.plot,
         xlab=xlab,
         ylab=ylab,
         panel.first=grid(),
         col=density.cols,
         log="xy",
         pch=compa.pch,
         cex=plot.cex,
         ...)
    
    ## Mark the alpha thresholds
    if (plot.colors) {
      abline(v=alpha, lty="solid", col="blue", lwd=1)
      abline(h=alpha, lty="solid", col="#00BB00", lwd=1)
    } else {
      abline(v=alpha, lty="dashed", col="black", lwd=1)
      abline(h=alpha, lty="dashed", col="black", lwd=1)      
    }
    
    # Add a legend
    if (plot.colors) {
      legend.col=c("white", "white", "blue", "#00BB00", "red", "black", "white", "white")
    } else {
      legend.col="black"
    }

    legend(legend.corner,
           c(paste("Total:", nrow(frame.to.plot)),
             paste(score.name, "<=", alpha),
             paste(xlab, ":", result$nb.signif.score1),
             paste(ylab, ":", result$nb.signif.score2),
             paste("both:", result$nb.signif.both),
             paste("none:", result$nb.signif.none),
             paste("Jaccard:", signif(digits=2, result$jaccard))
           ),
           bty="o",bg="white",
           pch=c(-1,-1,3,4,8,1,-1),
           cex=legend.cex,
           col=legend.col)
    
  }
}
