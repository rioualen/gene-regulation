################################################################
## A set of functions for plotting ROC curves, and computing the area
## under the curve (AUC)

################################################################
## Plot an empty frame for the ROC curve:
## - axes from 0 to 1
## - diagonals
## - grid every 0.1
#' @title Plot an empty frame to draw a ROC curve
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Plot an empty frame to draw a ROC curve
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param main="ROC curve"  Main title for the plot
#' @param xlab="False Positive Rate (FPR)"  Label for the abcsissa
#' @param ylab="True Positive Rate (TPR)" Label for the ordinate
#' @param xlim=c(0,1)  Limits for the abcsissa
#' @param ylim=c(0,1)  Limits for the ordinates
#' @param grid.bars=seq(from=0,to=1,by=0.1) Grid bars
#' @param plot.limits=TRUE Plot a black rectangle with the limits.
#' @param plot.diagonal=TRUE Plot the diagonal to indicate random expectation.
#' @param ... Additional parameters are passed to plot().
#' @export
plotROCframe <- function(main="ROC curve",
                         xlab="FPR",
                         ylab="Sensitivity",
                         xlim=c(0,1),
                         ylim=c(0,1),
                         grid.bars = seq(from=0,to=1,by=0.1),
                         plot.limits=TRUE,
                         plot.diagonal=TRUE,
                         ...) {

  ## Plot the frame
  plot(0, 
       0, 
       xlab=xlab, 
       ylab=ylab, 
       xlim=xlim, 
       ylim=ylim, 
       type="n",
       main=main, ...)

  ## Plot grid every 10%
  abline(v=grid.bars,col="#DDDDDD")
  abline(h=grid.bars,col="#DDDDDD")

  ## Plot limits
  if (plot.limits) {
    abline(h=0:1,col="#888888")
    abline(v=0:1,col="#888888")
  }
  
  ## Highlight the diagonal
  if (plot.diagonal) {
    abline(a=0,b=1,col="#888888")
    #  abline(a=1,b=-1,col="#BBBBBB")
  }
}


#' @title ROC statistics
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute the statistics required to draw a ROC curve, 
#' + some additional validation statistics
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param neg.scores
#' @param pos.scores
#' @param neg.total=length(neg.scores)
#' @param pos.total=length(pos.scores)
#' @param sort.decreasing=F  Indicate whether the score should be sorted by
#' decreasing order (this option should be used when lower scores are 
#' considered to be more significant, e.g. P-value, E-value, FDR).
#' @export
ROCstatistics <- function(neg.scores,
                          pos.scores,
                          neg.total=length(neg.scores),
                          pos.total=length(pos.scores),
                          sort.decreasing=F) {
  all.scores <- sort(unique(c(neg.scores, pos.scores)),decreasing=sort.decreasing)
  n.scores <- length(all.scores)
  
  ROC.table <- data.frame(row.names=all.scores,
                          score=all.scores,
                          pos.counts=rep(0,n.scores)
  )
  
  ## Distribution of positive scores
  pos.counts <- table(pos.scores)
  ROC.table$pos.counts <- 0
  ROC.table[names(pos.counts),"pos.counts"] <- pos.counts
  ROC.table$pos.counts.icum <-  rev(cumsum(rev(ROC.table$pos.counts)))
  ROC.table$pos.iCDF <-  ROC.table$pos.counts.icum/pos.total
  
  ## Distribution of negative scores
  neg.counts <- table(neg.scores)
  ROC.table$neg.counts <- 0
  ROC.table[names(neg.counts),"neg.counts"] <- neg.counts
  ROC.table$neg.counts.icum <-  rev(cumsum(rev(ROC.table$neg.counts)))
  ROC.table$neg.iCDF <-  ROC.table$neg.counts.icum/neg.total
  
  return(ROC.table)
}


#' @title Compute the area under the curve
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute the area under the curve from a table returned by ROCstatistics
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param ROC.table A result from ROC.table
#' @param ... Additional parameters are passed to XXXX().
#' @export
AUCfromROCtable <- function(ROC.table) {
  x <- c(0, sort(ROC.table$neg.iCDF), 1)
  y <- c(0, sort(ROC.table$pos.iCDF), 1)
  n <- length(x)
  
  auc <- t(as.matrix(abs(x[2:n] - x[1:n-1])))%*%as.matrix(abs((y[2:n] + y[1:n-1])/2))
  auc <- as.vector(auc)
  return(auc)
}


#' @title Plot a ROC curve
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Plot a ROC curve from a table returned by ROCstatistics. 
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param ROC.table, ## A data fame resulting from the function ROCstatistics()
#' @param main='ROC'
#' @param xlab='negative scores'
#' @param ylab='positive scores'
#' @param add=F If true, the ROC curve is added on a pre-existing plot
#' @param ... Additional parameters are passed to lines().
#' @export
plot.ROCstatistics <- function(ROC.table, ## A data fame resulting from the function ROCstatistics()
                               main='ROC',
                               xlab='negative scores',
                               ylab='positive scores',
                               add=F, ## If true, the ROC curve is added on a pre-existing plot
                               ...
) {
  if (!add) { plotROCframe(main=main,xlab=xlab,ylab=ylab) }
  x <- c(0, sort(ROC.table$neg.iCDF), 1)
  y <- c(0, sort(ROC.table$pos.iCDF), 1)
  lines(x,y,...)
}


################################################################
## Add labels on a ROC curve, for a given set of FPR values
labels.for.fpr <- function(ROC.table,
                           neg.values=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4,0.5)) {  
  label.frame <- data.frame()
  for (neg in neg.values) {
    i <- which(ROC.table$neg.iCDF < neg)
    i <- i[1]
    pos <- ROC.table[i,"pos.iCDF"]
    score <- ROC.table[i,"score"]
    label.frame <- rbind(label.frame,
                         data.frame(i,pos,neg,score))
  }
  segments(label.frame$neg, 0,label.frame$neg,label.frame$pos,col="darkred",lty="dotted")
  text(label.frame$neg,0,labels=round(label.frame$neg*100),pos=3,font=2,col="darkred",srt=90)
  segments(0,label.frame$pos,label.frame$neg,label.frame$pos,col="darkred",lty="dotted")
  text(0,label.frame$pos,labels=round(label.frame$pos*100),adj=c(1,0.5),font=2,col="darkred",bg='white',bty="n")
  text(label.frame$neg,label.frame$pos,labels=label.frame$score,adj=c(0,1),font=2,col="darkred")
  return(label.frame)
}


################################################################
## Old name maintained for backward compatiility
plot.ROC.curves <- function (...) {
  plot.ROC.curves.from.files(...)
}


## ##############################################################
## Plot a series of ROC curves from a set of data files.
## Each input file must be a tab-delimted text file with 3 columns:
## - score
## - frequency of true positive above the score (e.g.regulons)
## - frequencies of false positive above the score (e.g. random gene selections)
##
## One separate curve is plotted for each data file. This is
## convenient to compare performances of a program with different
## parameters, or to compare the performances obtained with different
## programs
plot.ROC.curves.from.files <- function(data.files, ## A list of files
                                       line.type = "l", ## line type for the plots
                                       main='ROC curve -  regulons versus random gene selections',
                                       xlab='random gene selections',
                                       ylab='regulons',
                                       ... ## Other parameters are passed to the plot() function
) {
  
  ## Plot the background lines
  plotROCframe(main=main,xlab=xlab,ylab=ylab,...)
  #   plot(c(0,1),
  #        c(0,1),
  #        type='l',
  #        panel.first=grid(col='#000000'),
  #        main=main,
  #        xlab=xlab,
  #        ylab=ylab,
  #        ...
  #        )
  #   lines(c(1,0),c(0,1),type='b',col='#BBBBBB')
  
  i <- 0
  data.colors <- vector()
  
  for (file.data in data.files) {
    ## Assign a color to each data file
    i <- i+1
    data.colors[file.data] <- i
    
    ## Read the data file
    setwd(dir.compa); x <- read.table(file.data,header=T)
    names(x) <- c('sig', 'regulons', 'random.genes')
    
    if (sum(x$random.genes) < sum(x$regulons)) {
      lines(x$random.genes, x$regulons,
            type=line.type,
            col=data.colors[file.data],
            lwd=2,
      )
    } else {
      lines(1- x$random.genes, 1- x$regulons,
            type=line.type,
            col=data.colors[file.data],
            lwd=2,
      )
    }
  }
  legend(0.5,0.5,legend=data.files,col=data.colors,lwd=2,bty="n")
  
}
