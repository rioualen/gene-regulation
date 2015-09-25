#' @title Descriptive statistics on each row of the input matrix
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute descriptive parameters (central tendency, dispersion)
#' for each row of a matrix or data frame.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#'
#' @return
#' A data.frame with one row per row of the input matrix, and one column
#' per computed statistics.
#' \describe{
#'  \item{mean}{mean}
#'  \item{median}{median}
#'  \item{sd}{standard deviation}
#'  \item{var}{variance}
#'  \item{iqr}{inter-quartile range}
#'  \item{sdiqr}{standardized IQR, which can serve as robust estimator of the standard deviation.}
#' }
#'
#' @examples
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#' data(denboer2009.group.labels)    ## Load phenotypic data
#'
#' stats.per.row <- statsPerRow(denboer2009.expr)
#' names(stats.per.row)
#' @export
statsPerRow <- function(x) {
  stats.per.row <- data.frame(
    mean = apply(x,1,mean),
    #     mean.goi = apply(x[,goi.columns],1,mean),
    #     mean.other = apply(x[,other.columns],1,mean),
    median = apply(x,1,median),
    #     median.goi = apply(x[,goi.columns],1,median),
    #     median.other = apply(x[,other.columns],1,median),
    sd = apply(x,1,sd),
    #     sd.goi = apply(x[,goi.columns],1,sd),
    #     sd.other = apply(x[,other.columns],1,sd),
    var = apply(x,1,var),
    #     var.goi = apply(x[,goi.columns],1,var),
    #     var.other = apply(x[,other.columns],1,var),
    iqr = apply(x,1,IQR)
    #     iqr.goi = apply(x[,goi.columns],1,IQR),
    #     iqr.other = apply(x[,other.columns],1,IQR),
  )
  stats.per.row$sdiqr = stats.per.row$iqr/(2*qnorm(3/4))
  #   stats.per.row$sdiqr.goi = stats.per.row$iqr.goi/(2*qnorm(3/4))
  #   stats.per.row$sdiqr.other = stats.per.row$iqr.other/(2*qnorm(3/4))
  row.names(stats.per.row) = row.names(x)

  return(stats.per.row)
}


