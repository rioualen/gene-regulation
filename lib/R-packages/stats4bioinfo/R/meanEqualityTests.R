#' @title Multiple test for mean equality
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Run different tests of mean equality (parametric or not)
#' on each row of a data table, and return a table summarizing the results
#' returned by each test.
#' 
#' The function also runs test to check the conditions of applicability 
#' for the mean equality tests: 
#' (1) Test of normality: Shapiro; 
#' (2) Tests of variance equality: F-test  (parametric), 
#'     Levene (non-parametric), 
#'     Brown-Forsythe (non-parametric, robust estimators).
#' 
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#' @param g       A vector describing group assigment
#' (length should equal the number of columns of the input matrix)
#' @param goi     Group of interest. Required for the two-groups tests (Student, Welch, Wilcoxon)
#' @param selected.tests=NA    Selection of a subset of tests to run.
#' @param alpha   Significance threshold (will be applied on corrected p-values)
#' @param verbosity level of verbosity (print messages during execition)
#'
#' @return
#' \describe{
#'   \item{param}{Parameters of the analysis}
#'   \item{selected.tests}{Tests ran. If no test was specified, all possible tests are selected.}
#'   \item{nb.per.group}{A sorted vector indicating the number of individuals (columns) per group.}
#'   \item{stats.per.row}{data.frame summarizing some descriptive statistics + the significance scores (p-value, e-value, fdr) returned by the
#'   different tests for each row of the input matrix.}
#' }
#'
#' @examples
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#' data(denboer2009.amp)    ## Load absent/marginal/present calls to filter input set
#'
#' ## Select subtypes represented by at least 30 samples
#' verbose("Selecting subtypes with at least 30 samples")
#' samples.per.subtype <- table(denboer2009.pheno$sample.labels)
#' selected.subtypes <- names(samples.per.subtype)[samples.per.subtype >= 30]
#' selected.samples <- denboer2009.pheno$sample.labels %in% selected.subtypes
#' 
#' ## Define group labels and group of interest
#' g <- denboer2009.pheno[selected.samples, "sample.labels"]
#' goi <- "Bh" ## Select one cancer subtype as group of interest
#'
#' ## Filter out genes called absent in most samples
#' verbose("Selecting probeset present in at least 30 samples")
#' selected.probesets <- apply(denboer2009.amp == "P", 1, sum) >= 30
#' x <- denboer2009.expr[selected.probesets, selected.samples]
#' verbose(paste("Matrix size:", nrow(x), "x", ncol(x)))
#' 
#' ## Run several mean equality tests on each row of the data set.
#' ## We restrict it to the first probesets for the demo.
#' diff.results <- meanEqualityTests(x=x, g=g, goi=goi, 
#'     selected.tests=c("welch", "wilcoxon"), verbosity=1)
#'
#' ## Return the parameters of the analysis
#' print(diff.results$param)
#' print(diff.results$selected.tests)
#' print(diff.results$nb.per.group)
#'
#' ## Compare p-values returned by 2 tests
#' welch.vs.wilcoxon <- meanEqualityTests.compareTwoTests(
#'    diff.results, test1="welch", test2="wilcoxon")
#'
#' ## Display the names of the result field stat.per.row
#' names(diff.results$stats.per.row)
#'
#' ## Draw a color Volcano plot for Welch test results
#' meanEqualityTests.plotVolcano(diff.results, test="welch", legend.cex=0.7, plot.cex=0.5)
#'
#' ## Draw a grayscale Volcano plot for Welch test results
#' meanEqualityTests.plotVolcano(diff.results, test="welch", legend.cex=0.7, plot.cex=0.5, plot.colors=FALSE)
#' @export
meanEqualityTests <- function(x,
                              g,
                              goi = NULL,
                              selected.tests = NULL,
                              alpha = 0.05,
                              verbosity=0) {
  
  
  ## Instantiate a list to return results
  result <- list()
  class(result) <- "meanEqualityTestResult"
  
  ## Ceck if selected tests were defined
  supported.tests.nogroup <- c("shapiro")
  supported.tests.2groups <- c("vartest", "levene2g", "brownforsythe2g", "student", "welch", "wilcoxon")
  supported.tests.ngroups <- c("bartlett", "levene", "brownforsythe", "anova", "kruskal")
  supported.tests <- c(supported.tests.nogroup,
                       supported.tests.2groups,
                       supported.tests.ngroups)
  if (is.null(selected.tests)) {
    selected.tests <-supported.tests
  }
  
  ## Check if selected tests are supported
  invalid.tests <- selected.tests[!selected.tests %in% supported.tests]
  if (length(invalid.tests) > 0) {
    print (paste("Invalid test", invalid.tests))
    stop("Selected tests contain invalid values.")
  }
  
  ## Check if the selected test require to specify a group of interest
  selected.tests.2groups <- selected.tests[selected.tests %in% supported.tests.2groups]
  if (length(selected.tests.2groups) > 0) {
    if (is.null(goi)) {
      print(paste(selected.tests.2groups, "requires to define a group of interest (goi)."))
      stop("You need to define a group of interest for some selected tests.")
    }
  }
  
  ## Compute the number of groups
  nb.groups <- length(unique(g))
  
  ## Indicate parameters in the result
  result$param <- c(ncol=ncol(x),
                    nrow=nrow(x),
                    nb.groups=nb.groups,
                    goi=goi,
                    alpha=alpha)
  result$selected.tests <- selected.tests

  ## Compute the number of samples per group
  result$nb.per.group <- as.vector(table(g))
  names(result$nb.per.group) <-  names(table(g))
  result$nb.per.group <- sort(result$nb.per.group, decreasing=TRUE)
  
  ## Compute row-wise statistics on the whole matrix
  verbose(paste(sep="", "Computing row-wise statistics (", 
                nrow(x), " features x ",
                ncol(x), " individuals)"), 1)
  
  stats.per.row <- statsPerRow(x)
  
  ## Define data structures required for 2 groups analysis (GOI versus other)
  if (!is.null(goi)) {
    
    ## Define the labels "goi versus other"
    verbose(paste("Defining", goi, "vs other labels"), 2)
    goi.vs.others.labels <- g
    goi.vs.others.labels[g != goi] <- "other"
    ##    verbose(table(goi.vs.others.labels), 1)
    
    
    ## Count the number of samples (columns of expression matrix) corresponding to the GOI
    goi.columns <- g==goi
    goi.nb <- sum(goi.columns)
    
    ## Count the number of samples (columns of expression matrix) corresponding to the other subtypes
    other.columns <- g!=goi
    other.nb <- sum(other.columns)
    
    ## Compute stats per row for the group of interest and for others.
    stats.per.row <- cbind(stats.per.row,
                           goi=statsPerRow(x[,goi.columns]),
                           other=statsPerRow(x[,other.columns]))
    names(stats.per.row)
    
    ## Check if the GOI exists in assigned classes
    if (sum(goi.vs.others.labels == goi) == 0) {
      stop(paste(sep="", "No instance of the group of interest (", goi, ") in the classes."))
    }
  }
  
  ###################################################
  ## Run a Shapiro-Wilk normality test on each
  ## row of the input matrix.
  if ("shapiro" %in% selected.tests) {
    verbose(paste(sep="", "Running Shapiro test of normality (", 
                  nrow(x), " features x ",
                  ncol(x), " individuals)"), 1)
    
    return.shapiro.p.value <- function(x){ shapiro.test(x)$p.value }
    stats.per.row$shapiro.p.value <- apply(x, 1, return.shapiro.p.value)
    
    ## Run a multiple testing correction on Shapiro p-values
    result$shapiro.multicor <- multipleTestingCorrections(p.value=stats.per.row$shapiro.p.value)
    stats.per.row$shapiro.fdr <- result$shapiro.multicor$multitest.table$fdr
    stats.per.row$shapiro.e.value <- result$shapiro.multicor$multitest.table$e.value
  }
  
  
  
  ###################################################
  ## Run a 2-group test of homoscedasticity (F variance test) on each
  ## row of the input matrix.
  if ("vartest" %in% selected.tests) {
    verbose(paste(sep="", "Running F variance test (",goi," vs other)"), 1)
    
    ## Run F variance test of homoscedasticity
    return.vartest.p.value <- function(x,g, ...){
      return(var.test(x=x[g==g[1]], y=x[g!=g[1]],...)$p.value)
    }
    
    stats.per.row$vartest.p.value <- apply(x,
                                        1,
                                        return.vartest.p.value,
                                        g=goi.vs.others.labels)
    
    result$vartest.multicor <- multipleTestingCorrections(p.value=stats.per.row$vartest.p.value)
    stats.per.row$vartest.fdr <- result$vartest.multicor$multitest.table$fdr
    stats.per.row$vartest.e.value <- result$vartest.multicor$multitest.table$e.value
  }

  ###################################################
  ## Run levene test of homoscedasticity on one row. 
  ## This function is used with various paramters to obtain: 
  ## - 2 groups variance equality test (by sending "GOI vs other" as group labels)
  ## - Brown-Forsythe robust equality test (by setting location to median rather than mean)
  return.levene.p.value <- function(x,g, ...){levene.test(x, g,...)$p.value}
  
  
  ###################################################
  ## Run Levene's robust test of homoscedasticity with 2 groups (GOI vs other) on each
  ## row of the input matrix.
  if ("levene2g" %in% selected.tests) {
    verbose(paste(sep="", "Running Levene 2-groups variance equality test (",goi," vs other)"), 1)
    library(lawstat)
    
    stats.per.row$levene2g.p.value <- apply(x,
                                       1,
                                       return.levene.p.value,
                                       g=goi.vs.others.labels,
                                       location="mean")
    
    result$levene2g.multicor <- multipleTestingCorrections(p.value=stats.per.row$levene2g.p.value)
    stats.per.row$levene2g.fdr <- result$levene2g.multicor$multitest.table$fdr
    stats.per.row$levene2g.e.value <- result$levene2g.multicor$multitest.table$e.value
  }
  
  
  ###################################################
  ## Run Brown-Forsythe's robust test of homoscedasticity with 2 groups (GOI vs other) on each
  ## row of the input matrix.
  if ("brownforsythe2g" %in% selected.tests) {
    verbose(paste(sep="", "Running Brown-Forsythe 2-groups variance equality test (",goi," vs other)"), 1)
    library(lawstat)
    
    stats.per.row$brownforsythe2g.p.value <- apply(x,
                                         1,
                                         return.levene.p.value,
                                         g=goi.vs.others.labels,
                                         location="median")
    
    result$brownforsythe2g.multicor <- multipleTestingCorrections(p.value=stats.per.row$brownforsythe2g.p.value)
    stats.per.row$brownforsythe2g.fdr <- result$brownforsythe2g.multicor$multitest.table$fdr
    stats.per.row$brownforsythe2g.e.value <- result$brownforsythe2g.multicor$multitest.table$e.value
  }
  
  
  ###################################################
  ## Run Levene's multi-group test of homoscedasticity on each
  ## row of the input matrix. Levene's test reputed more robust 
  ## to non-normality than Bartlett's test.
  if ("levene" %in% selected.tests) {
    verbose(paste(sep="", "Running Levene variance equality test (", nb.groups ," groups)"), 1)
    
    
    stats.per.row$levene.p.value <- apply(x,
                                       1,
                                       return.levene.p.value,
                                       g=g,
                                       location="mean")
    
    result$levene.multicor <- multipleTestingCorrections(p.value=stats.per.row$levene.p.value)
    stats.per.row$levene.fdr <- result$levene.multicor$multitest.table$fdr
    stats.per.row$levene.e.value <- result$levene.multicor$multitest.table$e.value
  }
  
  
  ###################################################
  ## Run Brown-Forsythe's multi-group test of homoscedasticity on each
  ## row of the input matrix. Brown-Forsythe's test reputed more robust 
  ## to non-normality than Levene's test.
  if ("brownforsythe" %in% selected.tests) {
    verbose(paste(sep="", "Running Brown-Forsythe variance equality test (", nb.groups ," groups)"), 1)
    
    
    stats.per.row$brownforsythe.p.value <- apply(x,
                                       1,
                                       return.levene.p.value,
                                       g=g,
                                       location="median")
    
    result$brownforsythe.multicor <- multipleTestingCorrections(p.value=stats.per.row$brownforsythe.p.value)
    stats.per.row$brownforsythe.fdr <- result$brownforsythe.multicor$multitest.table$fdr
    stats.per.row$brownforsythe.e.value <- result$brownforsythe.multicor$multitest.table$e.value
  }
  
  
  
  
  ###################################################
  ## Run a Bartlett multi-group test of homoscedasticity on each
  ## row of the input matrix. Note that Bartlett test is very sensitive 
  ## to non-normality. 
  if ("bartlett" %in% selected.tests) {
    verbose(paste(sep="", "Running Bartlett variance equality test (", nb.groups ," groups)"), 1)
        
    ## Run bartlett test of homoscedasticity
    return.bartlett.p.value <- function(x,g, ...){bartlett.test(x, g,...)$p.value}

    stats.per.row$bartlett.p.value <- apply(x,
                                         1,
                                         return.bartlett.p.value,
                                         g=g)

    result$bartlett.multicor <- multipleTestingCorrections(p.value=stats.per.row$bartlett.p.value)
    stats.per.row$bartlett.fdr <- result$bartlett.multicor$multitest.table$fdr
    stats.per.row$bartlett.e.value <- result$bartlett.multicor$multitest.table$e.value
  }




  ## Define a function to apply t.test in parallel, which will be used
  ## for both Student and Welch tests
  return.t.test.p.value <- function(x,y, ...){t.test(x[y==goi],
                                                  x[y!=goi],
                                                  alternative="two.sided",
                                                  ...)$p.value}



  ###################################################
  ## Run Student t.test to each row
  if ("student" %in% selected.tests) {
    verbose(paste(sep="", "Running Student t-test of mean equality (",goi," vs other)"), 1)

    ## Run the t.test and store the result in a column of the stats table
    # student.p.value.column <- paste(sep="", "student.p.value.", goi, ".vs.others" )
    stats.per.row$student.p.value <- apply(x,
                                        1,
                                        return.t.test.p.value,
                                        y=goi.vs.others.labels,
                                        var.equal=TRUE)

    ## Compute the FDR for row-wise Student test
    result$student.multicor <- multipleTestingCorrections(p.value=stats.per.row$student.p.value)
    stats.per.row$student.fdr <- result$student.multicor$multitest.table$fdr
    stats.per.row$student.e.value <- result$student.multicor$multitest.table$e.value
  }


  ###################################################
  ## Run Welch t.test to each row
  if ("welch" %in% selected.tests) {
    verbose(paste(sep="", "Running Welch t-test of mean equality (",goi," vs other)"), 1)

    ## Run Welch t.test and store the result in a column of the stats table
    stats.per.row$welch.p.value <- apply(x,
                                      1,
                                      return.t.test.p.value,
                                      y=goi.vs.others.labels,
                                      var.equal=FALSE)

    ## Apply multiple testing corrections
    result$welch.multicor <- multipleTestingCorrections(p.value=stats.per.row$welch.p.value)
    stats.per.row$welch.fdr <- result$welch.multicor$multitest.table$fdr
    stats.per.row$welch.e.value <- result$welch.multicor$multitest.table$e.value

  }


  ###################################################
  ## Run Wilcoxon test to each row
  if ("wilcoxon" %in% selected.tests) {
    verbose(paste(sep="", "Running Wilcoxon non-parametric test of mean equality (",goi," vs other)"), 1)

    ## Define a function to apply wilcoxon test to each row of a table
    return.wilcoxon.p.value <- function(x,y, ...) {
      wilcox.test(x[y==goi],
                  x[y!=goi], conf.int=FALSE,
                  alternative="two.sided",
                  ...)$p.value
    }

    ## Run wilcoxon t.test and store the result in a column of the stats table
    stats.per.row$wilcoxon.p.value <- apply(x,
                                         1,
                                         return.wilcoxon.p.value,
                                         y=goi.vs.others.labels)

    ## Apply multiple testing corrections
    result$wilcoxon.multicor <- multipleTestingCorrections(p.value=stats.per.row$wilcoxon.p.value)
    stats.per.row$wilcoxon.fdr <- result$wilcoxon.multicor$multitest.table$fdr
    stats.per.row$wilcoxon.e.value <- result$wilcoxon.multicor$multitest.table$e.value

  }



  ###################################################
  ## Run ANOVA to each row
  if ("anova" %in% selected.tests) {
    verbose(paste(sep="", "Running ANOVA test of mean equality (", nb.groups ," groups)"), 1)


    ## Define a function to apply ANOVA test to each row of a table
    return.anova.p.value <- function(x,y, ...) {
      ## Build a linear model and run anova() on it
      return(anova(lm(x ~ y))[1,"Pr(>F)"])
    }

    ## Run ANOVA and store the result in a column of the stats table
    stats.per.row$anova.p.value <- apply(x,
                                      1,
                                      return.anova.p.value,
                                      y=g)

    ## Apply multiple testing corrections
    result$anova.multicor <- multipleTestingCorrections(p.value=stats.per.row$anova.p.value)
    stats.per.row$anova.fdr <- result$anova.multicor$multitest.table$fdr
    stats.per.row$anova.e.value <- result$anova.multicor$multitest.table$e.value
    ## mulitpleTestingCorrections.plotPvalDistrib(anova.multicor)
  }

  ###################################################
  ## Run Kruskal-Wallis non-parametric test to each row
  if ("kruskal" %in% selected.tests) {
    verbose(paste(sep="", "Running Kruskal-Wallis non-parametric test of mean equality (", nb.groups ," groups)"), 1)
    ## Define a function to apply Kruskal-Wallis test to each row of a table
    return.kruskal.p.value <- function(x,y, ...) {
      return(kruskal.test(expr ~group, data=data.frame(expr=x, group=y))$p.value)
    }


    ## Run Kruskall-Wallis test and store the result in a column of the stats table
    stats.per.row$kruskal.p.value <- apply(x,
                                        1,
                                        return.kruskal.p.value,
                                        y=g)

    ## Apply multiple testing corrections
    result$kruskal.multicor <- multipleTestingCorrections(p.value=stats.per.row$kruskal.p.value)
    stats.per.row$kruskal.fdr <- result$kruskal.multicor$multitest.table$fdr
    stats.per.row$kruskal.e.value <- result$kruskal.multicor$multitest.table$e.value

  }


  result$stats.per.row <- stats.per.row
  return(result)
}


# ## TO DO: ADAPT THE PRINT METHOD FOR THE RESULTS
# print.meanEqualityTestsResult <- function(mnEqlT.result) {
#   if (class(mnEqlT.result) == "meanEqualityTestResult") {
#     to.print <- cat("BOUM")
#   } else {    
#     warning(class(mnEqlT.result))
#     stop("Invalid class for print.meanEqualityTestsResult: should be a result of print.meanEqualityTests()")
#   }
# 
#   return(to.print)
# }


#' @title Return a contingency table for several tests with a meanEqualityTests result
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Compare the significant features returned by two or more tests in 
#' an object returned by meanEqualityTest(), and plot the result in the form of a Venn diagram.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param mnEqlT.result   A result of stats4bioinfo::meanEqualityTests()
#' @param selected.tests=NULL  Tests to compare. If NULL, all available tests are compared.
#' @param signif.score="fdr"    Score on which the alpha threshold should be applied (Supported: fdr, p.value, e.value)
#' @param alpha=0.01           Significance threshold (will be applied on corrected p-values)
#' @param plot.venn=FALSE       Plot a Venn diagram with the result.
#' @param plot.tree=FALSE       Plot a hierarchical tree indicating relationships between tests.
#' @param hclust.method="average"  Hierarchical clustering method (passed to hclust).
#' @param plot.heatmap=FALSE  Plot a heatmap (in gray scale) highlighting the similarity between feature lists.
#' @param ...  Additional arguments are passed to plot().
#'
#' @return
#' A list containing the following fields:
#' \describe{
#'   \item{selected.per.test}{A dataframe with Boolean values indicating the status 
#'   (significant or not) for each row (feature) of the input dataset.}
#'   \item{contingency}{A contingency table indicating the number of features (rows) 
#'   passing the significance threshold for each pair of tests. The diagonal of this 
#'   contingency table indicates the number of significant features for each individual 
#'   test.}
#' }
#'
#' @examples
#' ## We forward to meanEqualityTests() for examples of utilization
#' example(meanEqualityTests)
#'
#' ## Compare 2-group mean equality test results with 1% threshold on FDR
#' test.compa.2groups <- meanEqualityTests.compareSignifFeatures(diff.results,
#'                                                     selected.tests=c("student", 
#'                                                                      "welch", 
#'                                                                      "wilcoxon"),
#'                                                     signif.score="fdr",
#'                                                     alpha=0.01,
#'                                                     plot.venn=TRUE)
#' summary(test.compa.2groups)
#' print(test.compa.2groups$contingency)
#'
#' ## Compare multi-group mean equality test results with a threshold of 1 on the 
#' ## e-value (accept one false positive per analysis)
#' test.compa.ngroups <- meanEqualityTests.compareSignifFeatures(diff.results,
#'                                                     selected.tests=c("anova", 
#'                                                                      "kruskal"),
#'                                                     signif.score="e.value",
#'                                                     alpha=1)   
#'                                                     
#' ## Compare all the mean equality test results with a threshold of 1 on the 
#' ## e-value (accept one false positive per analysis)
#' test.compa <- meanEqualityTests.compareSignifFeatures(diff.results,
#'                                                       selected.tests=c("student", 
#'                                                                        "welch", 
#'                                                                        "wilcoxon",
#'                                                                        "anova", 
#'                                                                        "kruskal"),
#'                                                       signif.score="e.value",
#'                                                       alpha=1,
#'                                                       plot.venn=FALSE,
#'                                                       plot.tree=FALSE,
#'                                                       hclust.method="average",
#'                                                       main="Common significant probesets",
#'                                                       xlab="Differential method",
#'                                                       ylab="Jaccard distance",
#'                                                       plot.heatmap=TRUE)   
#'                                                       
#'                                                                                                                                                         )
#' summary(test.compa)
#' print(test.compa$contingency)
#' 
#' 
#' @export
meanEqualityTests.compareSignifFeatures <- function(mnEqlT.result,
                                                    selected.tests=NULL,
                                                    signif.score="fdr",
                                                    alpha=0.05,
                                                    plot.venn=FALSE,
                                                    plot.tree = FALSE,
                                                    hclust.method = "average",
                                                    plot.heatmap = FALSE,
                                                     ...) {
#  attach(mnEqlT.result$stats.per.row)
  
  result <- list()
  
  
  ## Compute the score suffix for column selection
  score.suffix <- paste(sep='', '.', signif.score)
  
  ## Indicate the parameters in the returned object for tractability
  result$param <- data.frame(alpha=alpha,
                    signif.score = signif.score,
                    n.features = nrow(mnEqlT.result$stats.per.row)
                    )
  
  ## Identify all possible score columns, given the chosen significance score
  result$available.score.columns <- names(mnEqlT.result$stats.per.row)[grep(pattern = paste(sep='', score.suffix, '$'), x = names(mnEqlT.result$stats.per.row))]
  result$available.tests <- sub(pattern = score.suffix, replacement = '', result$available.score.columns)
  
  ## If selected tests are not specified, use all available tests
  if (is.null(selected.tests)) {
    selected.tests <- result$available.tests
  }
  result$selected.tests <- selected.tests
  
  ## Identify the columns containing the specified sig score for the selected tests.
  stat.columns <- paste(sep=".", selected.tests, signif.score)  
  result$selected.per.test <- data.frame(mnEqlT.result$stats.per.row[,stat.columns] <= alpha)
  #   names(result$selected.per.test) <- paste(sep=".", 
  #                                            names(result$selected.per.test), 
  #                                            alpha)
  names(result$selected.per.test) <- sub(pattern=score.suffix, rep='', names(result$selected.per.test))
  
  
  result$contingency <- t(as.matrix(result$selected.per.test)) %*% as.matrix(result$selected.per.test)
  result$union <- result$param$n.features - t(as.matrix(!result$selected.per.test)) %*% as.matrix(!result$selected.per.test)

  ## Compute matrices of Jaccard similarities and distances
  result$jaccard.sim <- result$contingency
  result$jaccard.sim[result$union != 0] <- result$contingency / result$union
  result$jaccard.dist <- 1 - result$jaccard.sim
  
#  detach(mnEqlT.result$stats.per.row)

  
  ## Draw Venn diagram if required
  if (plot.venn) {
     library(gplots)
     venn(result$selected.per.test)
#     result$venneuler <- venneuler(result$selected.per.test)
#     plot(result$venneuler)
#    library(VennDiagram)
  }
   
  ## Build a hierarchical tree
  result$tree <- hclust(d=as.dist(result$jaccard.dist), method = hclust.method)
  if (plot.tree) {
    plot(result$tree, ...)
  }

  ## Plot heatmap
  if (plot.heatmap) {
#     heatmap(result$jaccard.dist, 
#             hclustfun = hclust, 
#             Rowv=NA,
#             Colv=as.dendrogram(result$tree),
#             scale="none",
#             col=gray.colors(n=256, start = 0, end = 1))
    
    
    notecol <- "black"
    notecol[result$jaccard.sim > 0.5] <- "white"
    heatmap.2(result$jaccard.sim, 
              cellnote = round(digits=2,result$jaccard.sim),  # same data set for cell labels
              notecol=notecol,      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(12,9),     # widens margins around plot
              col=gray.colors(n=255, start = 1, end = 0),       # use on color palette defined earlier 
              #              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="col",     # only draw a row dendrogram
              Rowv=as.dendrogram(result$tree),
              Colv=as.dendrogram(result$tree),
              revC,TRUE,
              scale="none")
    
  }

  return(result)
}
