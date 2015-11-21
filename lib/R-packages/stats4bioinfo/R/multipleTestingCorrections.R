#' @title Multiple corrections for multiple testing
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Apply various multiple testing corrections on a list of P-values, report
#' all of them in a table with one row per test and one column per correction,
#' and optionally draw illutrative figures.
#'
#' In particular, we apply the elegant strategy from Storey & Tibshirani (2003)
#' to estimate the respective proportions of true and alternative hypotheses
#' among the tests.
#'
#'
#' @details
#' First version: 2012-02-10
#' Last modification: 2015-02
#'
#' @param p.value             list of p-values
#' @param names=NULL       Names associated to the p.value list (will be used as row.names for to result table)
#' @param lambda=0.5       lambda parameter for estimating pi0 = m0/m
#' @param alpha=0.05   Threshold of significance (alpha).
#' @param main='Multitesting corrections'
#' Prefix for the main title of the plots
#' @param run.qvalue=FALSE
#' For the sake of comparison, compute FDR and generate graphs
#' using Storey's qvalue library
#' (\url{http://genomics.princeton.edu/storeylab/qvalue/}).
#' @param plots=FALSE   If true, generate illutrsative plots
#' @param plot.pch=c(p.value=2,
#'                   fdr=4,
#'                   qval.0=3,
#'                   e.value=1,
#'                   fwer=20)   point type (character) associated to each statistics for the plot
#' @param plot.col=c(p.value='black',
#'                   fdr='blue',
#'                   qval.0='darkviolet',
#'                   e.value='darkbrown',
#'                   fwer='orange')
#' color for the plot
#' @param plot.elements=c("p.value", "fdr", "qval.0","e.value", "fwer")
#' Elements to draw on the plots (can be any subet of the default list).
#' @return
#' A list comprizing the following elements:
#'  \item{lambda}{Value provided for the input parameter lambda.}
#'  \item{alpha}{Value provided for the input parameter threshold (alpha).}
#'  \item{m}{total number of tests (= number of elements in the input list of p-values). m = m0 + m1.}
#'  \item{m0.est}{Estimation of the number of cases under null hypothesis.}
#'  \item{m1.est}{Estimation of the number of cases that do not comply with the null hypothesis (m1.est = m - m0.est).}
#'  \item{pi0}{Estimated proportion of trully null cases among all tests. pi0 = m0.est / m. }
#'  \item{nb.signif}{Number of tests declared significant above the specified alpha threshold}
#'  \item{multitest.table}{
#'  A table with one row per test, and one column per statistics (the fields described below).}
#'  \item{p.value}{P-value = probability to observe an at least as significant result under the null hypothesis. This column contains a replica of the input vector of p-values.}
#'  \item{rank}{rank by increasing P-value}
#'  \item{e.value}{E-value = expected number of false positives}
#'  \item{fdr}{False Discovery Rate = expected proportion of truly null hypotheses among the cases declared positives. The FDR is estimated using Storey and Tibshirani procedure.}
#'  \item{fwer}{Family-Wise Error Rate = probability to observe at least one FP among all the tests, assuming all of them are under null hypothesis.}
#'  \item{qval.0}{FDR assuming that all the tests are under null hypothesis (Benjamini-Hochberg approach).}
#' @references
#' 1. Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. JOURNAL-ROYAL STATISTICAL SOCIETY SERIES B.
#' 2. Storey, J.D. and Tibshirani, R. (2003) Statistical significance for genomewide studies. Proc Natl Acad Sci USA, 100, 9440–9445.
#' @examples
#' ## Generate a vector of 10,000 fake p-values
#' ## with a given proportion of trully null (m0=9,000) and
#' ## non-null (m1=1,000) cases.
#' my.p.values <- c(runif(n=9000, min=0, max=1),
#'                  runif(n=1000, min=0, max=1e-2))
#'
#' multitest.result <- multipleTestingCorrections(my.p.values)
#' attributes(multitest.result)
#'
#' ## In principle, m0 should be ~9,000
#' print(paste("m0.est =", multitest.result$m0.est))
#'
#' ## In principle, m1 should be ~1,000
#' print(paste("m1.est =", multitest.result$m1.est))
#'
#' ## In principle, pi0 should be ~0.9
#' print(paste("pi0 =", multitest.result$pi0))
#'
#' ## Plot P-value distributions
#' mulitpleTestingCorrections.plotPvalDistrib(multitest.result)
#'
#' ## Compare the different multiple testing corrections (Y axis)
#' ## versus the nominal p-value (X axis).
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result)
#'
#' ## Plot the number of significant features for different
#' ## multiple testing corrections
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result)
#' 
#' @export
multipleTestingCorrections <- function(p.value,
                                       names=NULL,
                                       lambda=0.5,
                                       alpha=0.05,
                                       main='Multitesting corrections',
                                       file.prefix=NULL,
                                       run.qvalue=FALSE,
                                       plots=FALSE
                                       ) {

  m <- length(p.value) ## Total number of tests (m = m0 + m1)

  ## Create a table for storing the results
  multitest.table <- data.frame(p.value=p.value)
  if (!is.null(names)) {
    rownames(multitest.table) <- names
  } else if (!is.null(names(p.value))) {
    rownames(multitest.table) <- names(p.value)
  }
  multitest.table$rank <- rank(p.value)

  ## E-value
  multitest.table$e.value <- p.value * m

  ## FDR
  multitest.table$fdr <- p.adjust(p.value, method="fdr")

  ## Family-wise error rate (fwer)
  multitest.table$fwer <- pbinom(q=0, size=m, prob=multitest.table$p.value, lower.tail=FALSE)
  #  multitest.table$fwer <- 1 - ( 1 - multitest.table$p.value)^m

  ## Compute q-value with Storey's library
  ##
  ## Beware: Storey's qvalue() function assumes that p-values are sorted
  ## by increasing order. However, for multilpeTestingCorrections we don't
  ## want to impose this constraint, because we want to be able to annotate
  ## Any user-provided table without changing the order ot hte rows.
  ## We fix this by sorting p-values before sending them to quvalue(),
  ## and re-ordering them by rank of the original p-value vector before
  ## adding them to the multitesting.result table.
  if (run.qvalue) {
    library(qvalue)
    qobj <- qvalue::qvalue(sort(multitest.table$p.value), lambda=lambda)
    multitest.table$qval.Storey <-  qobj$qvalues[rank(p.value)]
    #     if (plots) {
    #       ## Draw plots of the function qvalue::qplot
    #       qvalue::qplot(qobj)
    #       if (!is.null(file.prefix)) {
    #         dev.copy2pdf(paste(file.prefix, "_qvalue_plots.pdf"),
    #                      width=8,height=6);
    #       }
    #     }
  }
  
  ## q-value_0 (q-value corresponding to the case where all the features
  ## would be trully null).
  multitest.table$qval.0 <- multitest.table$p.value * m/multitest.table$rank
  ## multitest.table$qval.0[multitest.table$rank] <- cummax(multitest.table$qval.0[multitest.table$rank]) ## TO BE CHECKED
  message("ATTENTION: I MUST ADD THE CUMMAX FOR qval.0")
  
  ## Estimate numbers of null and alternative hypotheses
  m0.est <- min(m, sum(p.value >= lambda) / (1-lambda)) ## m0 cannot be larger than m
  m1.est <- m - m0.est
  pi0 <- m0.est/m

  #  print(c(lambda, m0.est))


  ## Prepare the result
  result <- list()
  if (run.qvalue) {
    result$qobj <- qobj
  }
  result$multitest.table <- multitest.table
  result$lambda <- lambda
  result$alpha <- alpha
  result$m <- m
  result$m0.est <- m0.est
  result$m1.est <- m1.est
  result$pi0 <- pi0
  result$nb.signif <- c(
    p.value=sum(multitest.table$p.value <= alpha),
    bonferoni=sum(multitest.table$e.value <= alpha),
    e.value.le.1=sum(multitest.table$e.value <= 1),
    fwer=sum(multitest.table$fwer <= alpha),
    qval.0=sum(multitest.table$qval.0 <= alpha),
    fdr=sum(multitest.table$fdr <= alpha)
  )
  if (run.qvalue) {
    result$nb.signif["qval.Storey"] = sum(multitest.table$qval.Storey <= alpha)
  }


  ################################################################
  ## Draw plots
  if (plots) {
    mulitpleTestingCorrections.plotPvalDistrib(multitest.result)
#    mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result)
#    mulitpleTestingCorrections.plotSignifFeatures(multitest.result)
  }

  return(result)
}


#' @title Multiple corrections for multiple testing
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Generate a plot that reproduces Fig 1 from Storey and Tibshirani (2003),
#' with some additional details, in order to illustrate the estimation of
#' the parameter PI0 = m0/m1.
#' @param multitest.result the list returned by the function multipleTestingCorrections().
#' @param main='Multitesting corrections'  main title of the plot
#' @param plot.legend=TRUE Plot a legend with some indicative numbers (m0, m1, pi0).
#' @param legend.corner="topright" corner wher the legend has to be placed.
#' @param legend.cex=1   Font size for the legend.
#' @param breaks=seq(from=0,to=1, by=0.05)    Breaks for the histogram
#' @param draw.lambda="arrow" Indicate the level of the lambda parameter. Supported 
#' representations: "arrow", "line", "none". 
#' @param draw.m0.line=TRUE Draw an horizontal line indicating the estimated m0 
#' (number of trully null features) per bin.
#' @param draw.mean.line=FALSE Draw a dashed horizontal line indicating the mean 
#' number of features per bin. The difference between this line and the "m0 per bin" 
#' line reflects the importance of truly alternative features.
#' @param col="#BBBBFF" Histogram background color (passed to hist()).
#' @param overlay=NULL  Boolean vector marking features of a special category 
#' (e.g. truly null features). A colored histogram will be printed on top of the 
#' main histogram, to indicate the number of features belonging to this group.
#' @param overlay.col="#CCCCCC" Color for the overlay histogram.
#' @param mean.line.col="black" Color to draw the line indicating the mean number of feature per bin.
#' @param m0.line.col="black" Color to draw the line indicating the estimated m0 per bin.
#' @param ...   Additional parameters are passed to hist()
#' @export
#' @examples
#' ## To obtain the input list (multitest.result), run the examples of
#' ## stats4bioinfo::multipleTestingCorrections().
#'
#' example(multipleTestingCorrections)
#'
#' ## Plot the p-value distribution + landmarks
#' mulitpleTestingCorrections.plotPvalDistrib(multitest.result, draw.lambda="line")
#'
mulitpleTestingCorrections.plotPvalDistrib <- function(multitest.result,
                                                       main='P-value distribution',
                                                       plot.legend=TRUE,
                                                       legend.cex=1,
                                                       legend.corner="topright",
                                                       breaks=seq(from=0,to=1, by=0.05),
                                                       draw.lambda = "arrow",
                                                       draw.m0.line=TRUE,
                                                       draw.mean.line=TRUE,
                                                       overlay=NULL,
                                                       col="#FFEEDD",
                                                       overlay.col="#CCCCCC",
                                                       mean.line.col="darkred",
                                                       m0.line.col="blue",
                                                       ...
                                                       ) {

  ## Plot the histogram
  distrib <- hist(multitest.result$multitest.table$p.value,
                  breaks=breaks, plot=TRUE,
                  xlab="nominal p-value",
                  main=main, col=col,
                  ...)
  m.per.bin <- sum(distrib$counts)/length(distrib$counts)
  m0.per.bin <- multitest.result$m0.est / length(distrib$count)
  
  ## Plot the histogram of overlay values, if specified
  if (!is.null(overlay)) {
    distrib.overlay <-  hist(multitest.result$multitest.table$p.value[overlay],
                             breaks=breaks, add=TRUE,col=overlay.col)
  }
  
  ## Draw horizontal line to show the expected frequency per histogram bar,
  ## and the frequency based on m0 estimation
  if (draw.mean.line) {
    abline(h=m.per.bin, col=mean.line.col, lty='dashed', lwd=2)
  }
  if (draw.m0.line) {
    abline(h=m0.per.bin, col=m0.line.col, lty='solid', lwd=2)
  }

  ## Draw an arrow to highlight the lambda parameter
  if (draw.lambda == "arrow") {
    arrows(multitest.result$lambda, m.per.bin*1.5, 
           multitest.result$lambda, m.per.bin*1.1, 
           angle = 30, length=0.05, code = 2, col=m0.line.col, lwd = 2)
  } else if (draw.lambda == "line") {
    abline(v=multitest.result$lambda, lty="dotted", col=m0.line.col, lwd=2)
  }

  ## Add a legend
  if (plot.legend) {
    legend(legend.corner, bty='o', bg='white', border=NA,
           cex=legend.cex,
           legend=(c(paste("lambda =", multitest.result$lambda),
                     paste("m =", multitest.result$m),
                     paste("per bin =", m.per.bin),
                     paste("m0.est =", round(multitest.result$m0.est)),
                     paste("m1.est =", length(multitest.result$multitest.table$p.value) - round(multitest.result$m0.est)),
                     paste("pi0 =", round(multitest.result$pi0, digits=3)))
           ))
  }
}

#' @title Plot corrected p-values versus nominal p-value
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Plot the different multi-testing corrected statistics as a
#' function of the nominal P-value.
#' @param multitest.result the list returned by the function multipleTestingCorrections().
#' @param main='Multitesting corrections'  main title of the plot
#' @param alpha=0.05   Threshold of significance (alpha).
#' @param plot.pch=c(p.value=2,fdr=4,qval.0=3,e.value=1,fwer=20)  Specific characters to distinguish the plotted statistics.
#' @param plot.col=c(p.value='#000000',fdr='#888888',qval.0='#666666',e.value='#BBBBBB',fwer='#444444') Specific colors or gray levels to distinguish the plotted statistics.
#' @param plot.elements=c("p.value","e.value","fwer","fdr","qval.0") Selection of elements to display on the plot.
#' @param plot.legend=TRUE  Plot a legend indicating the number of features declared significant with the alpha threshold on the selected statistics.
#' @param legend.corner="topleft" corner wher the legend has to be placed.
#' @param legend.cex=1   Font size for the legend.
#' @param xlab="p-value"
#' @param ylab="Multi-testing corrected statistics"
#' @param ...        Additional parameters are passed to plot()
#' @export
#' @examples
#' ## To obtain the input list (multitest.result), run the example of
#' ## stats4bioinfo::multipleTestingCorrections().
#' example(multipleTestingCorrections)
#' 
#' ## Plot all the multiple testing corrections at once
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result)
#'
#' ## Compare e-value and FWER
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result, plot.elements=c("e.value","fwer"))
#'
#' ## Compare e-value and FDR. 
#' ## This plot highlights the non-linear relationship between FDR and p-value.
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result, plot.elements=c("e.value","fdr"))
#'
#' ## Compare Benjamini-Hochberg (qval.0) and Storey-Tibshirani (fdr) estimates of FDR.
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result, plot.elements=c("qval.0","fdr"))
mulitpleTestingCorrections.plotCorrectedVsPval <- function (multitest.result,
                                                            main="Multitesting corrections",
                                                            xlab="p-value",
                                                            ylab="Multi-testing corrected statistics",
                                                            alpha=0.05,
                                                            legend.corner="topleft",
                                                            legend.cex=1,
                                                            plot.legend=TRUE,
                                                            plot.pch=c(p.value=2,
                                                                       fdr=4,
                                                                       qval.0=3,
                                                                       e.value=1,
                                                                       fwer=20),
                                                            plot.col=c(p.value='#000000',
                                                                       fdr='#888888',
                                                                       qval.0='#666666',
                                                                       e.value='#BBBBBB',
                                                                       fwer='#444444'),
                                                            plot.elements=c("p.value",
                                                                            "fdr",
                                                                            "qval.0",
                                                                            "e.value",
                                                                            "fwer"),
                                                            ...
                                                            ) {
  multitest.table <- multitest.result$multitest.table
  m <- multitest.result$m
  m0 <- multitest.result$m0
  Pi0 <- multitest.result$Pi0
  p.value.min = min(multitest.table$p.value)
  plot(multitest.table$p.value,
       multitest.table$p.value,
       main=main,
       xlab=xlab,
       ylab=ylab,
       log='xy',
       panel.first=c(
         abline(h=10^c(-10:3),col="#CCCCCC"),
         abline(v=10^c(-10:0),col="#CCCCCC"),
         abline(h=1,col="#666666")
       ),
       col=plot.col["p.value"], pch=plot.pch["p.value"],
       xlim=c(p.value.min, 1),
       ylim=c(p.value.min, m), type="n",
       ...)

  ## Plot the lines and prepare the legend
  legend.text <- vector()
  for (element in plot.elements) {
    lines(multitest.table$p.value, multitest.table[,element], col=plot.col[element], pch=plot.pch[element], type='p')
    n.below.alpha <- sum(as.vector(as.matrix(multitest.table[,element])) <= alpha)
    legend.text <- append(
      legend.text, 
      paste(sep="", element, " <= ", alpha, ": ", n.below.alpha)
    )
  }
  
  ## Show the alpha level
  abline(h=alpha, col='black', lty='dashed')
    
  ## Plot the legend
  if (plot.legend) {
    legend(legend.corner, 
           legend=legend.text,
           bg='#FFFFFF', bty="o", 
           cex=legend.cex,
           col=plot.col[plot.elements],
           pch=plot.pch[plot.elements],
    )
  }
  
}

#' @title Plot number of significant features as a function of the significance threshold.
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Plot the number of significant features as a function of the control
#' criterion (nominal p-value, e-value, fdr, ...).
#' @param multitest.result the list returned by the function multipleTestingCorrections().
#' @param main="Significant features"  main title of the plot
#' @param xlab="P-value derived statistics"
#' @param ylab="Significant features"
#' @param alpha=0.05   Threshold of significance (alpha).
#' @param plot.legend=TRUE  Plot a legend indicating the number of features declared significant with the alpha threshold on the selected statistics.
#' @param legend.corner="topleft" corner wher the legend has to be placed.
#' @param legend.cex=1   Font size for the legend.
#' @param plot.pch=c(p.value=2,e.value=1,fwer=20,fdr=4,qval.0=3)  Specific characters to distinguish the plotted statistics.
#' @param plot.col=c(p.value='#000000',e.value='#BBBBBB',fwer='#444444',fdr='#888888',qval.0='#666666') Specific colors or gray levels to distinguish the plotted statistics.
#' @param plot.elements=c("p.value","fdr","qval.0","e.value","fwer") Selection of elements to display on the plot.
#' @param ...        Additional parameters are passed to plot()
#'  
#' @examples
#' ## To obtain the input list (multitest.result), see the documentatio of
#' example(multipleTestingCorrections)
#'
#' ## Plot all the multiple testing corrections at once
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result)
#'
#' ## Compare e-value and FWER
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("e.value","fwer"))
#'
#' ## Compare e-value and FDR
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("fdr","e.value"))
#'
#' ## Compare Benjamini-Hochberg (qval.0) and Storey-Tibshirani (fdr) estimates of FDR
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("fdr","qval.0"))
#' @export
mulitpleTestingCorrections.plotSignifFeatures <- function (multitest.result,
                                                           main="Significant features",
                                                           xlab="Significance threshold",
                                                           ylab="Significant features",
                                                           alpha=0.05,
                                                           plot.legend = TRUE,
                                                           legend.corner="topleft",
                                                           legend.cex=1,
                                                           xlim=NULL,
                                                           plot.pch=c(p.value=2,
                                                                      fdr=4,
                                                                      qval.0=3,
                                                                      e.value=1,
                                                                      fwer=20),
                                                           plot.col=c(p.value='#000000',
                                                                      fdr='#888888',
                                                                      qval.0='#666666',
                                                                      e.value='#BBBBBB',
                                                                      fwer='#444444'),
                                                           plot.elements=c("p.value",
                                                                           "fdr",
                                                                           "qval.0",
                                                                           "e.value",
                                                                           "fwer"),
                                                           ...
                                                           ) {
  
  ## Retrieve the multitesting result table
  multitest.table <- multitest.result$multitest.table

  ## Define X limits
  p.value.min <- min(multitest.table$p.value) ## Minimal p-value
  m <- multitest.result$m ## Number of tests
  if (is.null(xlim)) {
    xlim <- c(p.value.min, m)
  }
  
  ## Draw the plot frame
  plot(x=multitest.table$p.value,
       y=multitest.table$rank,
       main=main,
       xlab=xlab,
       ylab=ylab,
       log='xy',
       panel.first=c(
         grid(equilogs=F, lty='solid', col='#CCCCCC'),
         abline(h=10^c(-11:3), col='#CCCCCC'),
         abline(h=1,col="#666666")
       ),
       col=plot.col["p.value"], pch=plot.pch["p.value"],
       xlim=xlim, type="n",
       ylim=c(1,nrow(multitest.table)),
       ...)
  
  ## Plot the lines and prepare the legend
  legend.text <- vector() ## Collect legend text
  n.below.alpha.values <- vector() ## Collect number of values below the alpha threshold for each stat
  for (element in plot.elements) {
    lines(x=multitest.table[,element], 
          y=multitest.table$rank, 
          col=plot.col[element], pch=plot.pch[element], type='p')
    n.below.alpha <- sum(as.vector(as.matrix(multitest.table[,element])) <= alpha)   
    n.below.alpha.values <- append(n.below.alpha.values, n.below.alpha)
    legend.text <- append(
      legend.text, 
      paste(sep="", element, " <= ", alpha, ": ", n.below.alpha)
    )
    #    segments(n.below.alpha, alpha, n.below.alpha, p.value.min, lty="dashed")
  }
  
  ## Plot landmark lines
  abline(v=alpha, col='black', lty='dashed')
  abline(h=n.below.alpha.values, lty="dashed")
  
  
  ## Plot the legend
  if (plot.legend) {
    legend(legend.corner, 
           legend=legend.text,
           bg='#FFFFFF', bty="o", 
           cex=legend.cex,
           col=plot.col[plot.elements],
           pch=plot.pch[plot.elements],
    )
  }
}
