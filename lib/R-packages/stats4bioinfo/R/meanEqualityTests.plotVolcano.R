#' @title Volcano plot from meanEqualityTest object
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Draw a volcano plot for a selected test of a meanEqualityTest result.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param mnEqlT.result   A result of stats4bioinfo::meanEqualityTests()
#' @param test="welch"    Statistical test for which the Volcano plot has to be drawn.
#' @param alpha           Significance threshold (will be applied on corrected p-values)
#' @param legend.corner="bottomright"   Position where legend should be placed.
#' @param plot.cex     Point size for the volcano (passed to graphics::plot() function).
#' @param legend.cex=1    Font size factor for legend. This allows to choose separate font sizes to draw the plot itself and the legend.
#' @param plot.colors=TRUE    Use colors for the plot   
#' @param ...             Additional parameters are passed to plot()
#'
#'
#' @examples
#' ## We forward to meanEqualityTests() for examples of utilization
#' example(meanEqualityTests)
#'
#' ## Draw a Volcano plot for Welch test results
#' meanEqualityTests.plotVolcano(diff.results, test="welch")
#'
#' ## Draw a Volcano plot for ANOVA results
#' meanEqualityTests.plotVolcano(diff.results, test="anova")
#'
#' @export
meanEqualityTests.plotVolcano <- function(mnEqlT.result,
                                          test="welch",
                                          main=paste(test, "Volcano plot"),
                                          alpha=0.05,
                                          legend.corner="bottomright",
                                          plot.cex=0.7,
                                          legend.cex = 1,
                                          plot.colors=TRUE,                                          
                                          ...) {
  
  
  ## Collect data for the plot
  n <- nrow(mnEqlT.result$stats.per.row)
  e.value <- mnEqlT.result$stats.per.row[,paste(sep="", test, ".e.value")]
  fdr <- mnEqlT.result$stats.per.row[,paste(sep="", test, ".fdr")]
  volcano.frame <- data.frame(mnEqlT.result$stats.per.row$goi.mean - mnEqlT.result$stats.per.row$other.mean, 
                              -log10(e.value))
  
  ## Define colors by density
  density.cols <- densCols(volcano.frame, nbin=256,
                           col=colorRampPalette(c("#BBBBBB", "#000000")))
  
  ## Add colors to highlight significance thresholds
  if (plot.colors) {
    density.cols[(fdr < alpha) & (e.value >= 1)] <- "#CCDD00"
    density.cols[(e.value < 1) & (e.value >= alpha)] <- "orange"
    density.cols[e.value < alpha] <- "red"
    legend.col <- c("black", "#CCDD00", "orange", "red")
  } else {
    density.cols[(fdr < alpha) & (e.value >= 1)] <- "#BBBBBB"
    density.cols[(e.value < 1) & (e.value >= alpha)] <- "#666666"
    density.cols[e.value < alpha] <- "black"
    legend.col <- c("black", "#BBBBBB", "#666666", "black")
  }
  
  ## Define point characters to highlight consistencies and inconsistencies
  volcano.pch <- rep(1, nrow(volcano.frame))
  volcano.pch[(fdr < alpha) & (e.value >= 1)] <- 3
  volcano.pch[(e.value < 1) & (e.value >= alpha)] <- 4
  volcano.pch[e.value < alpha] <- 8
  
  #   ## Define plot characters accordingto significance
  #   ## Highlight features passing the threshold on FDR
  #   points(volcano.frame[(fdr < alpha) & (e.value >= 1),], col="#CCDD00", pch=3, cex=0.7)  ## +
  #   
  #   ## Highlight features passing the lenient threshold on e-value
  #   points(volcano.frame[(e.value < 1) & (e.value >= alpha),], col="orange", pch=4,cex=0.7) ## x
  #   
  #   ## Highlight features passing the conservative threshold on e-value
  #   points(volcano.frame[e.value < alpha,], col="red", pch=8, cex=0.7) ## *
  #   
  
  plot(volcano.frame,
       panel.first=grid(),
       main=main,
       xlab="Mean difference (GOI vs others)",
       ylab="Significance = -log10(e.value)",
       col=density.cols,
       pch=volcano.pch,
       cex=plot.cex,
       ...
  )
  
  
  ## Draw some landmarks
  abline(v=0)
  if (plot.colors) {
    abline(h=0, col="orange")
    abline(h=-log10(alpha), col="red")
  } else {
    abline(h=0, lty="dotted")
    abline(h=-log10(alpha), lty="dashed")
  }
  legend(legend.corner,
         legend=c(paste("N=", n),
                  paste("FDR<=", alpha, ":", sum(fdr < alpha)),
                  paste("e-value<=", 1, ":", sum(e.value < 1)),
                  paste("e-value<=", alpha, ":", sum(e.value < alpha))),
         bty="o",bg="white",
         col=legend.col,
         pch=c(1, 3, 4, 8),
         cex=legend.cex,
  )
}

