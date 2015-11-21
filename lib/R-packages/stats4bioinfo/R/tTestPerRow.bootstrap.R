#' @title Estimate robustness of tTestPerRow by running a bootstrap procedure.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Estimate robustness of tTestPerRow by running a bootstrap procedure.
#' 
#' The standard bootstrap standard consists in anlayzing random selection of samples
#' of the same size as the original dataset, but with replacement. This sampling
#' is performed separately on the columns of the input table corresponding to the two
#'  subsets defined by the sample labels, and the sampled data table is analyzed 
#'  with tTestPerRow(). Note that this procedure ensures resampled groups of the same 
#'  size as the original groups, but induces redundancy within each group, and thus 
#'  violates the basic assumption of independence between samples that underlies the 
#'  t-test.
#' 
#' Alternatively, this function can also be used to perform a sub-sampling (selecting 
#' smaller sets of columns per group) with or without replacement. Sub-sampling without 
#' replacement presents the advantage of preserving the independence between elements 
#' of each group (as far as the elements of the original groups were independent), but 
#' reduces sample sizes, therefore resulting in a loss of power.
#' 
#' A trateoff can be achieved by the jacknife procedure, which consists in sampling, 
#' without replacement, subsets of sizes  n1-1 and n2-2, where ni is the sample size for
#' the ith group. For each group, each element is thus discarded in one of the 
#' resampling iterations. Note that the jacknife procedure is poorly efficient to 
#' temperate the effect of outliers, since a given outlier will be discarded in only
#' one of the ni jacknife iterations, but kept in the ni-1 other iterations. 
#' The current version does not implement the jacknife subsampling yet, but for 
#' sufficiently large sample sizes, an approximation can be achieved by sub-sampling 
#' with subset.sizes = c(n1 -1, n2 -1).
#' 
#' @details
#' First version: 2015-03
#' Last modification: 2015-05
#' 
#' @param x data frame containing one column per object and one row per feature
#' @param cl vector describing class membership
#' @param iterations Number of iterations of the bootstrap procedure
#' @param subset.sizes=NULL Named vector indicating the number os columns to 
#' select for each class name. If NULL (default), sizes are set equal to the 
#' class frequencies in the class membership vector (which only makes sense 
#' for the standard bootstrap procedure with replacement).
#' @param replace=TRUE Method for sampling the columns per class. 
#' The default is with replacement (bootstrap procedure). Setting replace=FALSE only 
#' makes sense for a sub-sampling procedure, with a specified subset.sizes vector 
#' indicating smaller values than the original class sizes.
#' @param discard.dup=FALSE Discard duplicated entries from the bootstrapped samples, 
#' to achieve a non-redundant bootstrap. 
#' In principle bootstrap is with replacement, which means that the resampled data 
#' has the same size as the original sample, but since it is drawn with replacement 
#' each element can be drawn 0, 1, 2, 3 or more times. This however creates a strong 
#' problem of dependency between samples, which violates the basic assumptions underlying
#' the t-test. To circumvent this problem, the option discard.dup filters out all duplicates 
#' from the bootstrapped list. The result should give similar results as a subsampling 
#' without replacement, of a size ~66% of the original sample size. A difference is that the 
#' number of elements will vary between iterations of the bootstrap, since it will depend 
#' on the particular number of duplicates drawn at random.
#' @param alpha=0.05 Significance threshold, which will be applied to count the 
#' number of rows passing the test with aplha on p-value, e-value and FDR, 
#' respectively.
#' @param support.quantile=0.75 Minimal percent of support required to declare a feature positive.
#' The default is to select features supported in 75\% of bootstrap iterations.
#' @param ... Additional parameters are passed to tTestPerRow().
#' 
#' @examples
#' ## Generate a random set with two samples from distinct 
#' ## populations A and B, characterized by and m rows (features),
#' ## among which m1 are under H1 (mA != mB) and m0 under 
#' ## H0 (mA = mB).
#' m1 <- 200 ## Number of features under H_1
#' m0 <- 200 ## Number of features under H_0
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
#' ## Check row-wise means per group
#' breaks=seq(from=-1, to=1+effect.size, by=0.025)
#' par(mfrow=c(2,1))
#' hist(apply(rand.2grp[,sample.labels=="A"],1, mean), breaks=breaks, main="Sample means, features under H0", xlab="mean per row", ylab="Rows", col="grey")
#' hist(apply(rand.2grp[,sample.labels=="B"],1, mean), breaks=breaks, main="Sample means, features under H1", xlab="mean per row", ylab="Rows", col="orange")
#' par(mfrow=c(1,1))
#' 
#' ################################################################
#' ## Apply Student t-test with classical bootstrap
#' ## (draw samples of same size as original groups, with replacement).
#' ## Set the support quantile to 0.95, in order to select only the features 
#' ## passing the alpha threshold in at least 95% of the iterations.
#' student.bootstrap <- tTestPerRow.bootstrap(x = rand.2grp, cl=sample.labels, 
#'    iterations=bootstrap.iterations, var.equal=TRUE, support.quantile=0.75, m2.minus.m1=TRUE)
#'    
#' ## Plot comparisons between p-values obtained for the two first bootstrap iterations
#' # x11(width=10, height=10)
#' par(mfrow=c(2,2))
#' plotPvalCompa(student.bootstrap$p.value[,1:2], score.name="p.value", alpha=0.05, main="Control on p-value", legend.corner="topleft")
#' plotPvalCompa(student.bootstrap$fdr[,1:2], score.name="fdr", alpha=0.05, main="Control on FDR", legend.corner="topleft")
#' plotPvalCompa(student.bootstrap$e.value[,1:2], score.name="e.value", alpha=1, main="Permissive control on e-value (<= 1)", legend.corner="topleft")
#' plotPvalCompa(student.bootstrap$e.value[,1:2], score.name="e.value", alpha=0.05, main="Stringent control on e-value", legend.corner="topleft")
#' par(mfrow=c(1,1))
#' 
#' 
#' ################################################################
#' ## Plot two histograms showing the distributions of support values for features under H1 and H0, respectively
#' # x11(width=5, height=8)
#' par(mfrow=c(2,2))
#' hist(student.bootstrap$stats.per.row$fdr.support[row.means.g2 == row.means.g1], 
#'      main="Bootstrap support for features under H0",
#'      xlab="Support",
#'      ylab="Number of features", col="grey",
#'      breaks=0:bootstrap.iterations)
#' abline(v=student.bootstrap$support.threshold, col="red", lwd=2)
#' hist(student.bootstrap$stats.per.row$fdr.support[row.means.g2 != row.means.g1], 
#'      main="Bootstrap support for features under H1",
#'      xlab="Support",
#'      ylab="Number of features", col="orange",
#'      breaks=0:bootstrap.iterations)
#' abline(v=student.bootstrap$support.threshold, col="red", lwd=2)
#' tTestPerRow.bootstrap.hist(student.bootstrap, plot.dcdf=TRUE, lwd=2,
#'      col="#BBCCFF", legend.corner="topright")
#' par(mfrow=c(1,1))
#' 
#' ################################################################
#' ## Subsampling: apply Student t-test on small-size sample subsets, 
#' ## drew witout replacement.
#' student.subsampled.n10 <- tTestPerRow.bootstrap(x = rand.2grp, cl=sample.labels, m2.minus.m1=TRUE,
#'     subset.sizes=c("A"=10,"B"=10), replace=FALSE, iterations=100, var.equal=TRUE, 
#'     support.quantile=0.75)
#' 
#' student.subsampled.n20 <- tTestPerRow.bootstrap(x = rand.2grp, cl=sample.labels, m2.minus.m1=TRUE,
#'     subset.sizes=c("A"=20,"B"=20), replace=FALSE, iterations=100, var.equal=TRUE, 
#'     support.quantile=0.75)
#' 
#' ################################################################
#' ## "Greedy" sub-sampling: sub-sampling without replacement and 
#' ## maximal subset sizes (n1-1, n2-1).
#' student.subsampled.greedy <- tTestPerRow.bootstrap(x = rand.2grp, cl=sample.labels, 
#'     subset.sizes=sample.sizes-1, replace=FALSE, iterations=100, var.equal=TRUE, 
#'     support.quantile=0.75)
#'                
#'                                              
#' ################################################################
#' ## Draw Volcano plots for the full dataset, bootstrapped 
#' ## and sub-sampled data.
#' # x11(width=10, height=10)
#' xlim <- c(
#'     round(min(unlist(student.bootstrap$stats.per.row$effect.size.min)), digits=1), 
#'     round(max(unlist(student.bootstrap$stats.per.row$effect.size.max)), digits=2))
#' ylim <- c(0, 1.1*ceil(-log10(min(unlist(student.bootstrap$fdr)))))
#' par(mfrow=c(3,2))
#' ## Classical volcano plot
#' tTestPerRow.plotVolcano(student.bootstrap$full.set.test, legend.corner="topleft", ylim=ylim, xlim=xlim,
#'     main=paste(sep="", "Volcano plot: ", "n1=", sample.sizes[1], ", n2=", sample.sizes[2]))
#' ## Volcano plot with confidence intervals
#' tTestPerRow.plotVolcano(student.bootstrap$full.set.test, legend.corner="topleft", ylim=ylim, xlim=xlim, plot.ci=TRUE,
#'     main=paste(sep="", "Confidence volcano: ", "n1=", sample.sizes[1], ", n2=", sample.sizes[2]))
#' ## Bootstrap volcano
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.bootstrap, legend.corner="topleft", 
#'     plot.sig.boxes=TRUE, plot.effect.boxes=FALSE, plot.rectangles=FALSE, ylim=ylim, xlim=xlim,
#'     main=paste(sep="", "Bootstrap volcano: ", "n1=", sample.sizes[1], ", n2=", sample.sizes[2]))
#' ## Subsampling volcano with subset.sizes=c(10,10)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.subsampled.n10, legend.corner="topleft", 
#'     plot.sig.boxes=TRUE, plot.effect.boxes=FALSE, plot.rectangles=FALSE, ylim=ylim, xlim=xlim,
#'     main=paste(sep="", "Subsampling volcano: ", "n1=", 10, ", n2=", 10))
#' ## Subsampling volcano with subset.sizes=c(20,20)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.subsampled.n20, legend.corner="topleft", 
#'     plot.sig.boxes=TRUE, plot.effect.boxes=FALSE, plot.rectangles=FALSE, ylim=ylim, xlim=xlim,
#'     main=paste(sep="", "Subsampling volcano: ", "n1=", 20, ", n2=", 20))
#' ## Greedy subsampling  with subset.sizes=c(n1-1,n2-1)
#' tTestPerRow.bootstrap.VolcanoBoxPlot(student.subsampled.greedy, legend.corner="topleft", 
#'     plot.sig.boxes=TRUE, plot.effect.boxes=FALSE, plot.rectangles=FALSE, ylim=ylim, xlim=xlim,
#'     main=paste(sep="", "Greedy subsampling: ", "n1=", sample.sizes[2]-1, ", n2=", sample.sizes[2]-1))
#' par(mfrow=c(1,1))
#'
#' 
#' @export
tTestPerRow.bootstrap <- function (x,
                                   cl,
                                   iterations = 100,
                                   subset.sizes = NULL,
                                   replace=TRUE,
                                   discard.dup=FALSE,
                                   alpha = 0.05,
                                   support.quantile=0.75,
                                   m2.minus.m1=FALSE,
                                   ...) {
  
  ## Initialize result tables
  resampled.p.value <- data.frame(matrix(NA, nrow=nrow(x), ncol=iterations)); row.names(resampled.p.value) <- row.names(x)
  resampled.e.value <- data.frame(matrix(NA, nrow=nrow(x), ncol=iterations)); row.names(resampled.e.value) <- row.names(x)
  resampled.fdr <- data.frame(matrix(NA, nrow=nrow(x), ncol=iterations)); row.names(resampled.fdr) <- row.names(x)
  resampled.t.obs <- data.frame(matrix(NA, nrow=nrow(x), ncol=iterations)); row.names(resampled.t.obs) <- row.names(x)
  resampled.effect.size <- data.frame(matrix(NA, nrow=nrow(x), ncol=iterations)); row.names(resampled.effect.size) <- row.names(x)
  
  
  
  ## Bootstrap: iterate t-tests with resampled values
  i <- 1
  for (i in 1:iterations) { 
    verbose(paste(sep="", "Bootstrap iteration: ", i, "/", iterations), 1)
    resampled.classes <- sampleClasses(cl, subset.sizes=subset.sizes, replace=replace, discard.dup=discard.dup)
    # hist(table(resampled.classes), breaks=0:(max(table(resampled.classes)))) ## Number of repetitions per sample
    # table(names(resampled.classes))
    # table(cl)
    sample.ttpr <- tTestPerRow(x = x[,resampled.classes], cl = names(resampled.classes), m2.minus.m1=m2.minus.m1) #, ...)
    resampled.p.value[,i] <- sample.ttpr$table$p.value
    resampled.e.value[,i] <- sample.ttpr$table$e.value
    resampled.fdr[,i] <- sample.ttpr$table$fdr
    resampled.t.obs[,i] <- sample.ttpr$table$t.obs
    resampled.effect.size[,i] <- sample.ttpr$table$means.diff
  }
  names(resampled.p.value) <- paste(sep="","p.value.", 1:iterations)
  names(resampled.e.value) <- paste(sep="","e.value.", 1:iterations)
  names(resampled.fdr) <- paste(sep="","fdr.", 1:iterations)
  names(resampled.t.obs) <- paste(sep="","t.obs.", 1:iterations)
  row.names(resampled.p.value) <- row.names(x)
  row.names(resampled.e.value) <- row.names(x)
  row.names(resampled.fdr) <- row.names(x)
  row.names(resampled.t.obs) <- row.names(x)
  
  
  ## Compute resampled.stats.per.row.per.row
  verbose(paste(sep="", "Computing stats per row from t-test results with resampled data"), 1)
  
  resampled.stats.per.row.per.row <- data.frame(
    "t.obs.mean" = apply(resampled.t.obs,1,mean),
    "t.obs.var" = apply(resampled.t.obs,1,var),

    "effect.size.mean" = apply(resampled.effect.size,1,mean),
    "effect.size.min" = apply(resampled.effect.size,1,min),
    "effect.size.q1" = apply(resampled.effect.size,1,quantile, 0.25),
    "effect.size.median" = apply(resampled.effect.size,1,median),
    "effect.size.q3" = apply(resampled.effect.size,1,quantile, 0.75),
    "effect.size.max" = apply(resampled.effect.size,1,max),
    "effect.size.var" = apply(resampled.effect.size,1,var),

    "sig.mean" = apply(-log10(resampled.e.value),1,mean),
    "sig.var" = apply(-log10(resampled.e.value),1,var),

    "p.value.min" = apply(resampled.p.value,1,min),
    "p.value.q1" = apply(resampled.p.value,1,quantile, 0.25),
    "p.value.median" = apply(resampled.p.value,1,median),
    "p.value.q3" = apply(resampled.p.value,1,quantile, 0.75),
    "p.value.max" = apply(resampled.p.value,1,max),
    "p.value.support" = apply(resampled.p.value <= alpha,1,sum),
    "p.value.support.quantile" = apply(resampled.p.value,1,quantile, support.quantile),
    
    "e.value.min" = apply(resampled.e.value,1,min),
    "e.value.q1" = apply(resampled.e.value,1,quantile, 0.25),
    "e.value.median" = apply(resampled.e.value,1,median),
    "e.value.q3" = apply(resampled.e.value,1,quantile, 0.75),
    "e.value.max" = apply(resampled.e.value,1,max),
    "e.value.support" = apply(resampled.e.value <= alpha,1,sum),
    "e.value.support.quantile" = apply(resampled.e.value,1,quantile, support.quantile),
    "e.value.le.1" = apply(resampled.e.value <= 1,1,sum),

    "fdr.min" = apply(resampled.fdr,1,min),
    "fdr.q1" = apply(resampled.fdr,1,quantile, 0.25),
    "fdr.median" = apply(resampled.fdr,1,median),
    "fdr.q3" = apply(resampled.fdr,1,quantile, 0.75),
    "fdr.max" = apply(resampled.fdr,1,max),
    "fdr.support" = apply(resampled.fdr <= alpha,1,sum),
    "fdr.support.quantile" = apply(resampled.fdr,1,quantile, support.quantile)
  )

  ## Run tTestPerRow on the full dataset
  verbose(paste(sep="", "Running t-test on the full dataset"), 1)
  full.ttpr <- tTestPerRow(x = x, cl = cl, m2.minus.m1=m2.minus.m1) #, ...)
  
  ## Collect results
  result<- list()
  result$cl <- cl
  result$iterations <- iterations
  result$alpha <- alpha
  result$samples.per.group <- table(cl)
  result$subset.sizes <- subset.sizes
  result$replace <- replace
  result$discard.dup <- discard.dup
  result$m2.minus.m1 <- m2.minus.m1
  result$support.quantile <- support.quantile
  result$support.threshold <- ceil(iterations*support.quantile) 
  result$nrow <- nrow(x)
  result$nsample <- ncol(x)
  result$effect.size <- resampled.effect.size
  result$p.value <- resampled.p.value
  result$e.value <- resampled.e.value
  result$fdr <- resampled.fdr
  result$t.obs <- resampled.t.obs
  result$stats.per.row <- resampled.stats.per.row.per.row
  result$full.set.test <- full.ttpr
  
  return(result)
}

#' @title Resample the elements of a named vector (typically class labels).
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Resample the elements of a named vector (typically class labels).
#' 
#' @param cl vector of class labels
#' @param subset.sizes=NULL named vector indicating the resampling size (values) for each class (names). If NULL, labels are resampled with the same frequency as in the input vector. 
#' @param replace=TRUE
#' 
#' @details
#' First version: 2015-04
#' Last modification: 2015-09
#' 
sampleClasses <- function(cl,
                          subset.sizes=NULL,
                          replace=TRUE,
                          discard.dup=FALSE) {
  
  class.names <- unique(cl) ## Attention, we need the names in the same order as in the vector of classes to avoid problems with the option m2.minus.m1

  ## Determine subset sizes
  if (is.null(subset.sizes)) {
    ## If no subset size is specified, use sample sizes as resampling sizes 
    ## (supposed to be used with replace=TRUE, to follow the classical bootstrap procedure).
    subset.table <- table(cl) ## Compute the frequency per class label
    subset.sizes <- as.vector(subset.table[class.names])
    names(subset.sizes) <- class.names
  } else {
    ## Check that subset sizes have names corresponding to sample labels
    if ((is.null(names(subset.sizes)))
        | (sum(!(names(subset.sizes) %in% unique(cl))) > 0)) {
      stop("Function sampleClasses(): subset.sizes must be vector, with entry names corresponding to class labels")
    }
    class.names <- names(subset.sizes)
  }
  
  cl.resampled <- vector()
  g <- 1 ## Instantiate group number for testing
  for (g in 1:length(subset.sizes)) {
    new.cl.selection <- sample(x=which(cl==class.names[g]), 
                               size=subset.sizes[g],
                               replace=replace)

    ## Discard duplicates for non-redundant bootstrap.
    if (discard.dup) {
      new.cl.selection <- unique(new.cl.selection)
    }
    names(new.cl.selection) <- rep(class.names[g],length(new.cl.selection))
    cl.resampled  <- 
      append(cl.resampled,new.cl.selection)
  }
  
  return(cl.resampled)
}
  
#' @title Plot an histogram with a bootstrap test result.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Plot an histogram from the result of tTestPerRow.bootstrap(). 
#' The histogram indicates the number of features (ordinate) with a given support
#' (abcissa), i.e. the number of features declared significant in exactly 
#' X bootstrap iterations. 
#' 
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' 
#' @param bootstrap.result    An object returned by the function tTestPerRow.bootstrap().
#' @param support.statistics="fdr"  Statistics for which the support has to be displayed.
#' Supported values: "fdr", "p.value", "e.value".
#' @param support.quantile=0.75 Minimal percent of support required to declare a feature positive.
#' The default is to select features supported in 75\% of bootstrap iterations.
#' @param xlab="Support"
#' @param ylab="Number of features"
#' @param main=paste(sep="", "Bootstrap support (", support.statistics,"<=", bootstrap.result$alpha, ")")
#' @param plot.dcdf=FALSE  Plot the decreasing CDF over the histogram, in "step" mode. 
#' This curve indicates the number of features supported by at least X bootstraps.
#' @param plot.legend = TRUE  Only valid if plot.dcdf is TRUE. Plot a legend indicating 
#' the number of features supported by a given number of the iterations. 
#' @param legend.corner="topright"  Corner to plo the legend.
#' @param lwd=1  Line width for the dcdf and support threshold lines.
#' @param ... Additional parameters are passed to hist().
#' 
#' @examples
#' ## First run the examples of tTestPerRow.bootstrap(), in order to get some result to plot
#' example("tTestPerRow.bootstrap")
#' 
#' ## Plot the histogram
#' tTestPerRow.bootstrap.hist(student.bootstrap, col="#BBBBBB")
#' 
#' ## Plot the histogram + the dcdf. 
#' tTestPerRow.bootstrap.hist(student.bootstrap, col="#BBBBBB", 
#'                            plot.dcdf=TRUE, lwd=2)
#'
#' ## Plot the histogram + the dcdf for the e-value. 
#' tTestPerRow.bootstrap.hist(student.bootstrap, support.statistics="e.value",
#'                            col="#BBBBBB", 
#'                            plot.dcdf=TRUE, lwd=2)
#' 
#' @export
tTestPerRow.bootstrap.hist <- function(bootstrap.result, 
                                       xlab = "Support",
                                       ylab = "Number of features",
                                       support.statistics = "fdr",
                                       support.quantile=bootstrap.result$support.quantile,
                                       main = paste(sep="", "Bootstrap support (", support.statistics,"<=", bootstrap.result$alpha, ")"), 
                                       plot.dcdf=FALSE, 
                                       plot.legend = TRUE,
                                       legend.corner="topright",
                                       lwd=1,
                                       ylim=NULL,
                                       ...) {
  support.values <- as.vector(as.matrix(bootstrap.result$stats.per.row[paste(sep=".", support.statistics, "support")]))
  
  ## Compute dcdf if required, and adapt Y axis limits
  if (plot.dcdf) {
    if (is.null(ylim)) {
      ylim <- c(0, nrow(bootstrap.result$stats.per.row))
    }
  }
  
  h <- hist(support.values, 
            breaks=0:bootstrap.result$iterations,
            main=main,
            xlab=xlab, 
            ylab=ylab, 
            ylim=ylim, ...)
  
  if (plot.dcdf) {
    dcdf <- rev(cumsum(rev(h$counts)))
    min.support <- floor(bootstrap.result$iterations * support.quantile)
    supported.features <- dcdf[min.support]
    lines(h$mids-1/2,dcdf, type="s", col="blue", lwd=lwd)
    arrows(min.support, 0, min.support, supported.features, length=0, col="brown", lwd=lwd)
    arrows(min.support, supported.features, 0, supported.features, length=0, col="brown", lwd=lwd)
    if (plot.legend) {
      legend(legend.corner, 
             col=c("blue", "brown"),
             legend=c(
               "dcdf",
               paste(sep="", "N(support >= ",min.support, ") = ", sum(support.values>=min.support))
             ), 
             lwd=lwd)
    }
  }
}


