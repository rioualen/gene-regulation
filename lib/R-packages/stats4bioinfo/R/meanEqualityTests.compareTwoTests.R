#' @title Compare the feature lists called significant by two tests from meanEqualityTest object.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Compare the feature lists  called significant by two tests 
#' (e.g. Welch versus Wilcoxon) from an object returned by meanEqualityTest(), 
#' and plot the result.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param mnEqlT.result   A result of stats4bioinfo::meanEqualityTests()
#' @param test1           First test
#' @param test2           Second test
#' @param signif.score    Score on which the alpha threshold should be applied (Supported: fdr, p.value, e.value)
#' @param alpha           Significance threshold (will be applied on corrected p-values)
#' @param ...             Additional parameters are passed to plotPvalCompa()
#'
#' @examples
#' ## We forward to meanEqualityTests() for examples of utilization
#' example(meanEqualityTests)
#'
#' ## Compare e-values between Welch and Wilcoxon tests
#' welch.vs.wilcoxon.e.value <- meanEqualityTests.compareTwoTests(diff.results,
#'                                                test1="welch",
#'                                                test2="wilcoxon",
#'                                                signif.score="e.value",
#'                                                plot.cex=0.5,legend.cex=0.8)
#'
#' ## Compare FDR between Welch and Wilcoxon tests, draw a grayscal plot
#' welch.vs.wilcoxon.fdr <- meanEqualityTests.compareTwoTests(diff.results,
#'                                                test1="welch",
#'                                                test2="wilcoxon",
#'                                                alpha=0.01,
#'                                                signif.score="fdr",
#'                                                plot.colors=FALSE,
#'                                                plot.cex=0.5,legend.cex=0.8)
#'
#' @export
meanEqualityTests.compareTwoTests <- function(mnEqlT.result,
                                              test1,
                                              test2,
                                              signif.score="fdr",
                                              alpha=0.05,
                                              ...) {
  result <- list()
  result$param <- c(test1=test1,
                    test2=test2,
                    alpha=alpha,
                    signif.score=signif.score
  )
  
  ## Collect data for the plot
  result$score.table <- data.frame(
    test1.score = mnEqlT.result$stats.per.row[,paste(sep=".", test1, signif.score)],
    test2.score = mnEqlT.result$stats.per.row[,paste(sep=".", test2, signif.score)]
  )
  
  ## Gather the data to plot
  frame.to.plot <- result$score.table[,c("test1.score", "test2.score")]
  
  result$compa <- plotPvalCompa(frame.to.plot, alpha=alpha, ...)
  
  return(result)
}



