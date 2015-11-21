#' @title T-test (Student or Welch) on each row of an input table.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Apply a t-test (Welch or Student) to each row of a multivariate
#' table. 
#' 
#' The estimators of central tendency and dispersion can either be computed 
#' with the moments (mean, standard deviation) or with robust esimators
#' (median, IQR).
#'
#' @details
#' First version: 2003-09
#' Last modification: 2015-02
#'
#' @param x                 A matrix or data frame
#' @param cl                A vector describing class assignment (length should equal the number of columns of the data table)
#' @param P.threshold=NA       p-value threshold. If specified, the result table only contains rows passing this threshold.
#' @param E.threshold=NA       e-value threshold. If specified, the result table only contains rows passing this threshold.
#' @param FDR.threshold=NA     Threshold on the False Discovery Rate (FDR). If specified, the result table only contains rows passing this threshold.
#' @param robust.est=F      Use robust estimators for central tendency and dispersion
#' @param verbosity=1       Level of verbosity
#' @param volcanoPlot=T     Draw a volcano plot.
#' @param alternative="two-sided"    Alternative hypothesis for the Welch test. Supported: "two.sided" (default), "less", "greater".
#' @param var.equal=F       Assume or nt var.equaliable. Apply Student test if true, Welch test if false.
#' @param m2.minus.m1=FALSE Compute the difference as (d=m2-m1) rather than (d=m1-m2). 
#' By default we set d=m1-m2 for the sake of consistency with stats::t.test().
#' @param group.col.name   Include group labels in the column name of the output table
#' @param alpha=0.05  Alpha value, used to define the confidence interval width
#' @param ...  Additional parameters are passed to the function tTestPerRow.plotVolcano()
#' 
#' @return
#' A data.frame with one row per test result, and one column per statistics.
#' @examples
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#' data(denboer2009.group.labels)    ## Load phenotypic data
#'
#' ## Print cancer types and associated group labels
#' print(data.frame(denboer2009.group.labels))
#'
#' ## Compute the number of samples per subtype of cancer (ALL)
#' sort(table(denboer2009.pheno$sample.labels), decreasing=TRUE)
#'
#' ## Create a vector with group labels per sample,
#' ## For the Welch test we compare one group of interest (e.g. Bh)
#' ## to all the other ones (labeled as "other").
#' goi <- "Bh" ## Group of interest
#' sample.groups <- denboer2009.pheno$sample.labels
#' sample.groups[sample.groups != goi] <- "other"
#'
#' ## Check number of samples per group
#' sort(table(sample.groups))
#'
#' ## Run Welch test on each row of the DenBoer dataset
#' welch.table <- tTestPerRow(x=denboer2009.expr, cl=sample.groups, var.equal=FALSE)
#' 
#' ## Draw a volcano plot with Welch result table
#' tTestPerRow.plotVolcano(welch.table, 
#'      main=paste(sep="", "Welch test : Den Boer (2009), ", goi, " vs others"),
#'      legend.corner = "topright")
#' 
#' ## Run Student test on each row of the DenBoer dataset
#' student.table <- tTestPerRow(x=denboer2009.expr, cl=sample.groups, var.equal=TRUE)
#'
#' ## Draw a volcano plot with Student result table
#' tTestPerRow.plotVolcano(student.table, 
#'             main=paste(sep="", "Student test : Den Boer (2009), ", goi, " vs others"),
#'             legend.corner="topleft")
#' 
#' ## Draw a volcano plot with confidence interval bars for Student result table
#' tTestPerRow.plotVolcano(student.table, 
#'             main=paste(sep="", "Student test : Den Boer (2009), ", goi, " vs others"),
#'             legend.corner="topleft", plot.ci=TRUE, alpha=0.05)
#' 
#' ## Compare e-values from Student and Welch tests
#' plotPvalCompa(data.frame(
#'    "Student"=student.table$e.value,
#'    "Welch"=welch.table$e.value), score="e-value", alpha=0.05)
#'    
#' ## Compare FDR from Student and Welch tests
#' plotPvalCompa(data.frame(
#'    "Student"=student.table$fdr,
#'    "Welch"=welch.table$fdr), score="FDR", alpha=0.05)
#'    
#' ## Confusion table between Student and Welch tests
#' table(student.table$fdr < 0.05, welch.table$fdr < 0.05) ## Lenient threshold on FDR
#' table(student.table$e.value < 0.05, welch.table$e.value < 0.05) ## Lenient threshold on E-value
#' table(student.table$e.value < 1, welch.table$e.value < 1) ## Intermediate threshold on E-value
#' 
#' @export
tTestPerRow <- function (x,
                         cl,
                         P.threshold=NA,
                         E.threshold=NA,
                         FDR.threshold=NA,
                         robust.est=F, ## use robust estimators for central tendency and dispersion
                         verbosity=1, ## Level of verbosity
                         volcanoPlot = FALSE, ## Draw a volcano plot
                         alternative = "two.sided", ## Supported: "two.sided", "less", "greater"
                         var.equal = FALSE, ## Student or Welch test
                         group.col.names=FALSE, ## Include the group label in the column names
                         m2.minus.m1=FALSE, ## Compute difference as m2 - m1 rather than m1 - m2
                         alpha=0.05, ## Alpha value, used to define the confidence interval width
                         ... ## additional parameters are passed to the function tTestPerRow.plotVolcano()
) {
  
  
  ## Dimensions of the data set
  n <- nrow(x)
  p <- ncol(x)
  
  ## check the dimensions of the class vector
  if (length(cl) != p) {
    stop (paste(sep='', 
                'Invalid column number (',p,')',
                'should equal the length of the class vector (',length(cl),')'))
  }

  ## Compute classical or robust estimators of the central tendency (mean) 
  ## and dispersion (variance sd).
  group.sizes <- table(cl)
  c <- unique(cl)
  if (length(unique(cl)) != 2) {
    stop("tTestPerRow: invalid vector of groups. Should contain exactly two different groups.")
  }

  ## Report starting time and parameters
  if (var.equal) {
    test.name <- "Student"
  } else {
    test.name <- "Welch"
  }
  verbose(paste(sep="", "tTestPerRow: ",test.name," t-test on ", nrow(x), " rows. ", 
                "Group sizes: ", group.sizes[1], ",", group.sizes[2], "."), 1)
  
  
  if (robust.est) {
    ## Actually robust estimators are problematic, I should suppress this option or send a strong warning !!!
    m1.est <- apply(x[,cl==c[1]],1,median,na.rm=T)
    m2.est <- apply(x[,cl==c[2]],1,median,na.rm=T)
    m.est <- apply(x,1,median,na.rm=T)    
    iqr1 <- apply (x[,cl==c[1]], 1, IQR, na.rm=T)
    iqr2 <- apply (x[,cl==c[2]], 1, IQR, na.rm=T)
    iqr <- apply (x, 1, IQR, na.rm=T)
    sd1.est <- iqr1/(qnorm(0.75) - qnorm(0.25))
    sd2.est <- iqr2/(qnorm(0.75) - qnorm(0.25))
    sd.est <- iqr/(qnorm(0.75) - qnorm(0.25))    
    var1.est <- sd1.est^2
    var2.est <- sd2.est^2
    var.est <- sd.est^2
    
  } else {
    m1.est <- apply(x[,cl==c[1]],1,mean,na.rm=T)
    m2.est <- apply(x[,cl==c[2]],1,mean,na.rm=T)
    m.est <- apply(x,1,mean,na.rm=T)    
    var1.est <- apply(x[,cl==c[1]],1,var,na.rm=T)
    var2.est <- apply(x[,cl==c[2]],1,var,na.rm=T)
    var.est <- apply(x,1,var,na.rm=T)
    sd1.est <- sqrt(var1.est)
    sd2.est <- sqrt(var2.est)
    sd.est <- sqrt(var.est)
  }
  
  n1 <- sum(cl == c[1])
  n2 <- sum(cl == c[2])
  
  ## 2012-05-02: restored means.diff = m1.est - m2.est for the sake
  ## of consistency with R function t.test()
  if (m2.minus.m1) {
    means.diff <- m2.est - m1.est
  } else {
    means.diff <- m1.est - m2.est
  }
  ##  hist(means.diff,breaks=50)

  
  ## Pooled variance and effect size
  var.pooled <- ((n1-1)*var1.est + (n2-1)*var2.est) / (n1+n2-2)
  sd.pooled <- sqrt(var.pooled)
  effect.size <- means.diff / sd.pooled
    
  ## Estimate the standard error of the difference and compute the degrees of freedom. 
  ## These two computations differ between Welch and Student tests.
  if (var.equal) {   
    ## Student test (assumption of equal variance)
    diff.stder <- (sd.pooled * (sqrt(1/n1+1/n2)))
    
    ## Calculate degrees of freedom with Students formula
    df <- (n1 + n2 - 2)

  } else {
    ## Welch test (assumption of unequal variance)
  
    ## Estimate the standard error of the difference
    diff.stder <- sqrt(var1.est/n1 + var2.est/n2)

    ## Calculate degrees of freedom with Welch's formula
    df <- (var1.est/n1 + var2.est/n2)^2 / ((var1.est/n1)^2/(n1-1) + (var2.est/n2)^2/(n2-1))
    
  }
  
  ## Compute the t statistics
  t.obs <- means.diff/diff.stder
  
  ## Calculate p-value and e-value
  #  p.value.normal.approx <- 2*pnorm(abs(t.obs),lower.tail=F)
  if (alternative == "greater") {
    p.value <- pt(t.obs,df,lower.tail=FALSE)
  } else if (alternative == "less") {
    p.value <- pt(t.obs,df,lower.tail=TRUE)
  } else if (alternative == "two.sided") {
    p.value <- 2*pt(abs(t.obs),df,lower.tail=F)
  } else {
    stop('Invalid alternative option for tTestPerRow(). Supported: "two.sided", "less", or "greater"')
  }
  e.value <- p.value*nrow(x)
  sig <- -log(e.value)/log(10)

  multi.corr <- multipleTestingCorrections(p.value, plots=FALSE)
  
  ## Compute the width of a confidence interval around the difference between means
  if (alternative == "two.sided") {
    critical.t <- qt(p=1-alpha/2, df=df)
  } else {
    critical.t <- qt(p=1-alpha, df=df)
  }
  ci.width <- 2*critical.t * diff.stder
  
  ## Build the result
  result.table <- data.frame(
    m1.est,
    m2.est,
    m.est,    
    means.diff,
    var1.est,
    var2.est,
    var.est, 
    var.pooled,
    sd1.est,
    sd2.est,
    sd.est,
    sd.pooled,
    effect.size,
    diff.stder,
    t.obs,
    df,
    p.value,
    e.value,
    sig)
  result.table$fwer <- multi.corr$multitest.table$fwer
  result.table$q.value <- multi.corr$multitest.table$qval.Storey
  result.table$fdr <- multi.corr$multitest.table$fdr
  result.table$rank <- multi.corr$multitest.table$rank
  result.table$ci.width <- ci.width
  #                       p.value.normal.approx,
  
  ## Adapt column headers to the type of estimators, and to the group names
  if (group.col.names) {
    if (robust.est) {
      names(result.table)[1:3] <- c(
        paste("median.", c[1], sep=""),
        paste("median.", c[2], sep=""),
        "medians.diff"
      )
    } else {
      names(result.table)[1:3] <- c(
        paste("mean.", c[1], sep=""),
        paste("mean.", c[2], sep=""),
        "means.diff"
      )
    }
    
    names(result.table)[4:7] <- c(
      paste("var.est.", c[1], sep=""),
      paste("var.est.", c[2], sep=""),
      paste("sd.est.", c[1], sep=""),
      paste("sd.est.", c[2], sep="")
    )
  }
  
  ## Filtering on p-value and e-value thresholds
  if (!is.na(P.threshold)) {
    result.table <- result.table[result.table$p.value < P.threshold,]
  }
  if (!is.na(E.threshold)) {
    result.table <- result.table[result.table$e.value < E.threshold,]
  }
  if (!is.na(FDR.threshold)) {
    result.table <- result.table[result.table$fdr < FDR.threshold,]
  }

  ## Report done time
  verbose("Multiple t-test done", 2)

  ## Draw a volcano plot
  if (volcanoPlot) {
    tTestPerRow.plotVolcano(result.table, ...)
  }
  
  ## Build the result object with result table + parameters of the analysis
  result <- list()
  result$table <- result.table
  result$alpha <- alpha
  result$p.threshold <- P.threshold
  result$e.threshold <- E.threshold
  result$fdr.threshold <- FDR.threshold
  result$alternative <- alternative
  result$var.equal <- var.equal
  result$m2.minus.m1 <- m2.minus.m1
  
  ##  plot(p.value.normal.approx,p.value,panel.first=grid(col='#0000ff'),log="xy")
  return (result)
}


#' @title Standard error plot
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Complement to the volcano plots generated by stats4bioinfo::tTestPerRow.plotVolcano()
#' Draw a plot that enables to compare
#' - standard error and significance
#' - mean difference, standard error and significance
#' @param multi.t Must be the data frame obtained with tTestPerRow()
#' @param col.points='#888888' Color(s) for the points
#' @param col.lines='blue'  Color for the line highlighting the significance threshold
#' @param col.grid=#CCCCCC   Grid color
#' @param col.positive='#00BB00' Color to highlight significant points (called positive)
#' @param ... Additional parameters are passed to the plot function
#' @examples
#' ## The object welch.table should have been generated with \code{\link{tTestPerRow}} (see examples in the doc)
#' stderPlot(welch.table)
#' @export
stderPlot <- function(multi.t,
                      col.points = '#888888',
                      col.grid = '#CCCCCC',
                      col.lines = 'blue',
                      col.positive='#00BB00',
                       ...
) {

  ## Draw standard error as a function of Welch significance
  plot(multi.t[,c('diff.stder','sig')],
       xlab='Standard error of mean difference',
       ylab='Welch significance (sig=-log10(e-value)',
       col=col.points,
       panel.first=grid(lty='solid',col=col.grid))
  abline(h=0,col=col.lines,lwd=2)
  points(multi.t[multi.t$sig > 0,c('diff.stder','sig')],col=col.positive)
}



stderPlotCircles <- function(multi.t,
                             col.points = '#888888',
                             col.grid = '#CCCCCC',
                             col.lines = 'blue',
                             col.positive='#00BB00',
                      ...  ## Plot standard error versus mean differences, with circles
                      ) {

  ## proportional to significance
  plot(multi.t[,c('means.diff', 'diff.stder')],
       cex=multi.t$sig-min(multi.t$sig),
       ylab='Standard error on the difference',
       xlab='Difference between the means',
       col=col.points,
       panel.first=grid(lty='solid',col=col.grid))
  abline(v=0, col=col.lines)
  points(multi.t[multi.t$sig > 0,c('means.diff', 'diff.stder')],
         cex=multi.t[multi.t$sig > 0,'sig']-min(multi.t$sig),
         col=col.positive)
}


# #' @title Compare two tTestPerRow results.
# #' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
# #' @description
# #' Compare two result tables produced by tTestPerRow
# #' @param ttpr.result1  First result produced by tTestPerRow
# #' @param ttpr.result2  Second result produced by tTestPerRow
# #' @param col.points='#888888' Color(s) for the points
# #' @param col.lines='blue'  Color for the line highlighting the significance threshold
# #' @param col.grid=#CCCCCC   Grid color
# #' @param col.positive='#00BB00' Color to highlight significant points (called positive)
# #' @param legend.corner="bottomright"   Position where legend should be placed.
# #' @param plot.cex     Point size for the volcano (passed to graphics::plot() function).
# #' @param ... Additional parameters are passed to the plot function
# #' @examples
# #' @export
# tTestPerRow.compare <- function(ttpr.result1,
#                                 ttpr.result2) {
#   
#   
# }
# 
# 
