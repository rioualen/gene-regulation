#' @title Generate a data frame where each row follows a custom random normal distribution.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description Generate a data frame where each row contains random normal values
#' with the same mean and standard deviation as the corresponding row in the 
#' input data frame.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#'
#' @return
#' A data frame of same dimensions as the input matrix/data frame, 
#' with random normal values generated in a row-wise way.
#'
#' @examples
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#'
#' ## Define group labels and group of interest
#' g <- denboer2009.pheno$sample.labels
#' goi <- "Bh" ## Select one cancer subtype as group of interest
#'
#' ## Generate row-wise normal values
#' denboer2009.rnorm <- rowwiseRandomNormal(x=denboer2009.expr)
#'
#' ## Run Welch test on the row-wise permuted values
#' rnorm.welch <- meanEqualityTests(denboer2009.rnorm, 
#'                                  g=denboer2009.pheno$sample.labels, 
#'                                  goi="Bh",
#'                                  selected.tests="welch"
#'                                  )
#'                
#' ## Draw volcano plot of Welch test result with the random values.
#' ## This should show more or less no significant features.
#' meanEqualityTests.plotVolcano(rnorm.welch, test="welch", main="Random normal values, Welch volcano")
#'
#' @export
rowwiseRandomNormal <- function(x) {
  verbose(paste(sep="", "Generating row-wise random normal matrix (", 
                nrow(x), " rows x ", ncol(x), " columns)"), 1)
  return.rnorm.one.row <- function(x){ rnorm(n=length(x), mean=mean(x), sd=sd(x)) }
  return(t(apply(x, 1, return.rnorm.one.row)))
}

#' @title Permute the values of each row of the input matrix or data frame.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description 
#' Random permutation of the values on each row of the input data frame.
#'
#' After permutation, each row thus contains exactly the same values as in 
#' the original expression matrix, but there should be no specific order or 
#' distinction between groups. 
#' 
#' Row-wise permuted matrices provide realistic negative controls for several 
#' approaches in microarray or RNA-seq analysis: differential expression, 
#' clustering, supervised classification.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#' @param ncol=ncol(x)   Number of columns of the sampled matrix  (passed as argument size to base::sample)
#' @param replace=FALSE   Sampling with replacement (passed to base::sample).
#'
#' @param ...     Additional parameters are passed to base::sample(). 
#' This extends the possible usages of stats4bioinfo::rowwiseSample().
#' For example, the argument replace=TRUE will perform a bootstrap of the input
#' table: each value of the input matrix can appear 0, 1 or more times in 
#' the corresponding row of the permuted matrix. 
#' 
#' @return
#' A data frame of same dimensions as the input matrix/data frame (unless the 
#' option size is used), with row-wise permuted values.
#'
#' @examples
#' 
#' ################################################################
#' ## Generate an artificial matrix with each row identical
#' x <- matrix(data=1:10, nrow=20, ncol=10, byrow=TRUE)
#' print(x)
#' 
#' ################################################################
#' ## Permute the values of each row of the matrix
#' print(rowwiseSample(x, replace=FALSE))
#' 
#' ################################################################
#' ## Sampling with replacement
#' print(rowwiseSample(x, replace=TRUE))
#' 
#' ################################################################
#' ## Sampling with a smaller number of columns
#' print(rowwiseSample(x, columns=4, replace=FALSE))
#' 
#' ################################################################
#' ## Sampling with a larger number of columns and replacement
#' print(rowwiseSample(x, columns=15, replace=TRUE))
#' 
#' ################################################################
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)     ## Load phenotype table
#'
#' ################################################################
#' ## Define group labels and group of interest
#' g <- denboer2009.pheno$sample.labels
#' group1 = "Bo" ## First cancer type for mean comparison
#' group2 = "Bt" ## Second cancer type for mean comparison
#' selected.samples <- (g == group1) | (g == group2)
#' selected.labels <- denboer2009.pheno$sample.labels[selected.samples]
#' print (paste("Groups", group1, "and", group2, ":", sum(selected.samples), "selected samples"))
#'
#' ## Run Welch test on the original values
#' denboer.welch <- meanEqualityTests(
#'     denboer2009.expr[selected.samples], 
#'     g=selected.labels, goi=group1,
#'     selected.tests="welch")
#'                  
#' ################################################################
#' ## Draw volcano plot of Welch test result on denboer2009 dataset
#' meanEqualityTests.plotVolcano(denboer.welch, test="welch", 
#'     main="Den Boer 2009, Welch volcano", 
#'     legend.corner="topright",)
#'     
#' ################################################################
#' ## Draw the distribution of p-values for Den Boer dataset. 
#' ## There is a strong peak at low p-values (<5%), which corresponds 
#' ## to the differentially expressed genes.
#' mulitpleTestingCorrections.plotPvalDistrib(
#'      denboer.welch$welch.multicor,
#'      main="Den Boer 2009, Welch p-values")
#'      
#' ################################################################
#' ## Permute the values of denboer2009 and run Welch test
#' permuted.profiles <- rowwiseSample(x=denboer2009.expr[selected.samples])
#'
#' perm.welch <- meanEqualityTests(permuted.profiles, 
#'    g=selected.labels, goi=group1,
#'    selected.tests="welch")
#'                  
#' ################################################################
#' ## Draw volcano plot of Welch test result with the permuted values
#' ## This should show more or less no significant features.
#' meanEqualityTests.plotVolcano(perm.welch, 
#'      test="welch", main="Permuted Den Boer 2009, Welch volcano",
#'      legend.corner="topright")
#'
#'                                            
#' ################################################################
#' ## Draw the distribution of p-values for the row-wise permuted matrix. 
#' ## This should give a flat distribution. However, the estimated proportion
#' ## of truly null (pi0) is lower than 1, because sometimes the permutation
#' ## creates unbalanced groups -> the randomized samples have different means. 
#' mulitpleTestingCorrections.plotPvalDistrib(
#'      perm.welch$welch.multicor, 
#'      main="Permuted Den Boer 2009, Welch p-values")
#' par(mfrow=c(1,1))
#'  
#' @export
rowwiseSample <- function(x,
                          columns=ncol(x),
                          replace=FALSE) {
  verbose(paste(sep="", "Generating row-wise permuted matrix (", 
                nrow(x), " rows x ", ncol(x), " columns)"), 1)
  
  x <- as.matrix(x)
  
  ## Define the function that will be applied to each row
  sampleOneRow <- function(row, 
                           columns, 
                           replace=FALSE) {
    permuted <- sample(as.vector(row), size=columns, replace=replace)
    return(permuted)
  }
  return(t(apply(x, 1, sampleOneRow, columns, replace)))
}

#' @title Permute rows in a group-balanced way
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description 
#' #' Random permutation of the values on each row of the input data frame, 
#' with a preservation of the proportions between 2-group labeled columns:
#' group of interest (GOI), and  "others", respecively.
#'
#' Each row thus contains exactly the same values as in the original expression
#' matrix, but there should be no specific distinction between groups.
#'
#' ATTENTION ! This procedure may give rise to surprizing bias. 
#' A priori it seemed to me that a balanced representation of the
#' original groups between the permuted samples would be a good 
#' idea, because I sometimes observed that permutation tests with small 
#' sample sizes would return too many significant results since sometimes 
#' the resampled groups contain different proportions of the original
#' groups. I thus implemented balanced permutation to suppress this 
#' effect. 
#' 
#' However, I noticed the opposite effect: when the effect size
#' is very strong in a given dataset, the balanced permuted set 
#' has an *under-representation* lof low p-values 
#' (e.g. 0 <= pval <= 30%), see example below. The cause of this 
#' surprizing behaviour is that the balanced permutations ensure
#' equality of the resampled group means, but if the original groups
#' have very different means, the resampled distributions are bimodal,
#' and have thus a high variance. The consequence is to artificially 
#' reduce the denominator of the t statistics ($t_{obs}$).
#'
#' 
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#' @param g       Group labels
#' @param goi     Group of interest. If not specified, the first label is take as group of interest.
#' In case g contains more than two distinct labels, these group labels are converted
#' to "GOI", and "other", respectively.
#'
#' @return
#' A data frame of same dimensions as the input matrix/data frame, with row-wise
#' group-balanced permuted values.
#'
#' @examples
#' ## Run example for rowwiseSample, in order to load 
#' ## the data and parameters
#' example(rowwiseSample)
#'
#' ## Permute the values of denboer2009
#' balanced.perm.profiles <- rowwisePermGroupBalanced(
#'     x=denboer2009.expr[selected.samples],
#'     g=selected.labels,
#'     goi=group1)
#'
#' ## Run Welch test on the row-wise permuted values
#' balanced.perm.welch <- meanEqualityTests(
#'     balanced.perm.profiles, 
#'     g=selected.labels, goi=group1,
#'     selected.tests="welch")
#'                  
#' ## Draw volcano plot of Welch test result with the permuted values, resp.
#' ## NOTE: we already see that the negative control is "too good": 
#' ## the highest significances are at -2 instead of 0.
#' meanEqualityTests.plotVolcano(balanced.perm.welch, 
#'    test="welch", 
#'    legend.corner='topright',
#'    main="Permuted Den Boer 2009, Welch volcano")
#' 
#' ## Plot p-value distribution for the balanced row-wise permuted dataset
#' ## NOTE: this plot clearly shows the bias of balanced permutation:
#' ## the low p-values (<= 30%) are under-represented because when the 
#' ## population means differe, the balanced resampling creates groups
#' ## with same expected mean, but increaseed variance.
#' mulitpleTestingCorrections.plotPvalDistrib(
#'    balanced.perm.welch$welch.multicor, legend.corner="bottomright",
#'    col='#FFBBBB')
#' 
#' @export
rowwisePermGroupBalanced <- function(x,
                                     g,
                                     goi=g[1]) {

  verbose(paste(sep="", "Generating group-balanced row-wise random normal matrix (", 
                nrow(x), " rows x ", ncol(x), " columns)"), 1)
  goi.vs.others.labels <- g
  goi.vs.others.labels[goi.vs.others.labels !=goi] = "other"

  ## Define the number of elements to sample in each group
  N <- ncol(x) ## Number of columns of input table
  n <- sum(g==goi)
  m <- N - n     ## Number of entries from other groups

  verbose(paste(sep="", "Group sizes: ", n, " versus ", m, " (N=",N, ")"), 1)
  
  ## Use a trick to randomly choose between ceiling and floor for the half of odd numbers
  n1.1 <- round((n + runif(1,-0.1,+0.1))^2/N) ## Number to samples from group 1 that remain in group 1
  n1.2 <- n - n1.1 ## Number of samples from group 2 assigned to group 1
  n2.1 <- n - n1.1 ## Number of samples from group 1 assigned to group 2
  n2.2 <- m - n1.2 ## Number of samples from group 2 remaining in group 2

  ## Define the function that will be applied to each row
  permuteOneRowGroupBalanced <- function(x, n) {
    ## We resample the sampled subgroups to break their order
    g1 <- sample(c(sample(1:n, size=n1.1), sample((n+1):N, size=n1.2)))
    g2 <- sample(c(sample(1:n, size=n2.1), sample((n+1):N, size=n2.2)))
    permuted <- x[c(g1, g2)]
    return(permuted)
  }

  return(t(apply(x, 1, permuteOneRowGroupBalanced, n=n)))
}


