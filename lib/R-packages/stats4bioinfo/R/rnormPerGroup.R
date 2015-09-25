#' @title Generate a data frame with random normal values sampled according to group-specific parameters.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description Generate a data frame with random normal values sampled according to 
#' group-specific parameters. Columns correspond to individuals belonging to different
#' groups characterized by specific means and/or standard deviations. Rows correspond 
#' to features.
#'
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#'
#' @param n       A vector indicating the number of columns per group.
#' @param mean    A vector indicating the mean per group. Must have the same length as n.
#' @param sd      A vector indicating the standard deviation per group. Must have the same length as n.
#' @param nrow    Number of rows (features) of the result data frame.
#' 
#' @return
#' A list with the following objects:
#' \describe{
#'   \item{x}{Data frame with the random numbers}
#'   \item{cl}{Vector with the class label of each column.}
#'   \item{mean.per.group}{Data frame with feature-wise means (rows) for each group (column).}
#'   \item{sd.per.group}{Data frame with feature-wise standard deviation (rows) for each group (column).}
#'   \item{exp.mean.per.col}{Vector with the expected means per column.}
#'   \item{mean.per.col}{Vector with the means per column in the result matrix.}
#'   \item{exp.sd.per.col}{Vector with the expected sds per column.}
#'   \item{sd.per.col}{Vector with the sds per column in the result matrix.}
#' }
#' A data frame of random normal values sampled with group-specific parameters.
#'
#' @examples
#' 
#' ################################################################
#' ## Small test: generate a matrix, composed of three groups of different sizes, means and sds.
#' small.rnorm <- rnormPerGroup(n=c(6,4,5), mean=c(-3,0,5), sd=c(2,1,4), nrow=10)
#' 
#' ## Check column means
#' plot(small.rnorm$exp.mean.per.col, small.rnorm$mean.per.col)
#' abline(a=0,b=1)
#' 
#' ## Check column sd
#' plot(small.rnorm$exp.sd.per.col, small.rnorm$sd.per.col)
#' abline(a=0,b=1)
#' 
#' ## Generate a wider matrix with
#' rnorm.result <- rnormPerGroup(n=c(100,100), mean=c(0, 0.5), sd=c(1,2), nrow=1000)
#' 
#' ## Check the means per column
#' boxplot(rnorm.result$mean.per.col ~ rnorm.result$cl)
#' 
#' ## Check the sd per column
#' boxplot(rnorm.result$sd.per.col ~ rnorm.result$cl, main="SD per group")
#' 
#' ## Check means per group
#' boxplot(rnorm.result$mean.per.group, main="Feature-wise mean per group")
#' 
#' ################################################################
#' ## Run Student test on each feature to check the power.
#' ## We chose equal sd to comply with the homoscedaticity assumption.
#' rnorm.result <- rnormPerGroup(n=c(50,50), mean=c(0, 0.5), sd=c(1,1), nrow=10000)
#' x.student <- tTestPerRow(x = rnorm.result$x, cl = rnorm.result$cl, var.equal=TRUE)
#' 
#' ## Plot histogram of the observed differences between groups
#' hist(x.student$table$means.diff, breaks=100, main="Effect size distribution", xlab="Effect size")
#' grid(lty="solid",col="#BBBBBB")
#' abline(v=0.5, col="blue", lwd=2)
#' 
#' ## Plot the histogram of p-values. 
#' ## Note: since all the data was generatd under H1, 
#' ## the distribution should be merely composed of low p-values.
#' hist(x.student$table$p.value, breaks=20, main="P-value distribution", xlab="p-value")
#' 
#' ## Plot the empirical power curve beta = f(alpha).
#' ## In this configuration where all features are under alternative hypothesis,
#' ## this corresponds to a Receiver-Operator Characterisitic (ROC) curve
#' ## With empirical TPR versus theoretical FPR.
#' plot(ecdf(x.student$table$p.value), 
#'    xlab=expression(FPR == alpha), ylab=expression(TPR == 1-beta),
#'    main=paste("Student ROC curve"), col="blue")
#' grid()
#' abline(v=c(0,1))
#' abline(h=c(0,1))
#' abline(a=0,b=1, lty="dashed")
#' 
#' 
#' @export
rnormPerGroup <- function(n, mean, sd, nrow) {
  
  ## Check consistency between  parameter vectors 
  group.nb <- length(n)
  verbose(paste("Generating random normal matrix with", group.nb, "groups"))
  if (length(mean) != group.nb) {
    stop(paste("rnormPerGroup() error: parameter 'mean' must be a vector of same length as parameter 'n'."))
  }
  if (length(sd) != group.nb) {
    stop(paste("rnormPerGroup() error: parameter 'sd' must be a vector of same length as parameter 'n'."))
  }
  column.mean <- rep(mean, n)
  column.sd <- rep(sd, n)
  ncol <- length(column.mean)
  group.names <- paste(sep="", "g", 1:group.nb)
  cl <- rep(group.names, n) ## Class names
  
  
  ## Generate the random values
  x <- data.frame(matrix(rnorm(n=nrow*ncol, mean = column.mean, sd = column.sd), byrow = TRUE, ncol=ncol))
  
  ## Set column names indicating the group
  column.names <- vector()
  mean.per.group <- data.frame(matrix(nrow=nrow, ncol=group.nb))
  sd.per.group <- data.frame(matrix(nrow=nrow, ncol=group.nb))
  for (g in 1:group.nb) {
    group.name <- group.names[g]
    column.names <- append(column.names, paste(sep=".", group.name, 1:n[g]))
    mean.per.group[,g] <- apply(x[,cl==group.name], 1, mean)
    sd.per.group[,g] <- apply(x[,cl==group.name], 1, sd)
  }
  colnames(x) <- column.names
  colnames(mean.per.group) <- group.names
  colnames(sd.per.group) <- group.names
  
  ## Generate the result object
  result <- list()
  result$x <- x
  result$cl <- cl
  result$mean.per.group <- mean.per.group  
  result$sd.per.group <- sd.per.group  
  result$exp.mean.per.col <- rep(mean,n)
  result$mean.per.col <- apply(x,2,mean)
  result$exp.sd.per.col <- rep(sd,n)
  result$sd.per.col <- apply(x,2,sd)
  
  return(result)
}  
  
