#' @title Add some columns to a functional enrichment table in order to ensure that we always have the same statistics.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Add some columns to a functional enrichment table in order to ensure that we always have the same statistics.
#'
#' @details
#' First version: 2015-08
#' Last modification: 2015-08.
#'
#' @param enrich.table Enrichment table from GOstats, clusterProfiler or similar programs.
#' @param pvalue.column Name or number of the column containing the p-values (this depends on the program).
#' @param thresholds=c("evalue"=1)
#' @examples
#'
#' @export
completeEnrichTable <- function(enrich.table,
                                pvalue.column = "Pvalue",
                                thresholds=c("evalue"=1, "qvalue"=0.05),
                                nb.tests=nrow(enrich.table),
                                select.positives=FALSE,
                                plot.adjust = TRUE
) {

  ## Ensure homogeneous column names exist in all results even if it duplicates the p-value column.
  if (pvalue.column != "pvalue") {
    enrich.table$pvalue <- enrich.table[,pvalue.column]
  }

  # Initialize the "positive" flag. It will be update for each selection criterion.
  enrich.table$positive <- 1
  # View(enrich.table)

  ## Compute the E-value
  enrich.table$evalue <- enrich.table[,pvalue.column] * nb.tests

  ## Sort enrichment table by increasing p-values
  enrich.table$pvalue.rank <- rank(enrich.table[,pvalue.column])
  enrich.table <- enrich.table[order(enrich.table[,pvalue.column], decreasing=FALSE),]

  ## q-value according to Benjamini-Hochberg method, which does not consider the prior proportions of null and non-null hypotheses
  enrich.table$qvalue <- enrich.table[,pvalue.column] * nb.tests/enrich.table$pvalue.rank
  # TO DO: cummax(x) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ## Compute the q-value according to Benjamini-Hochberg method

  ## Tag significant genes according to zero, one or more cutoff values
  upper.thresholds <- c("evalue"=TRUE, "qvalue"=TRUE, "pvalue"=TRUE, "intersect"=FALSE)
  for (s in names(thresholds)) {
    if (is.numeric(thresholds[s])) {
      selection.column <- paste(sep="_", s, thresholds[s])
      if (upper.thresholds[s]) {
        enrich.table[, selection.column] <- 1*(enrich.table[,s] < thresholds[s])
      } else {
        enrich.table[, selection.column] <- 1*(enrich.table[,s] > thresholds[s])
      }
      sum(enrich.table[, selection.column])
      verbose(paste(sep="", "\t",
                    sum(enrich.table[selection.column]), "/", nb.tests, " tests passed ",
                    s, " threshold=", thresholds[s]), 2)
      enrich.table$positive <- enrich.table$positive * enrich.table[,selection.column]
      verbose(paste(sep="", "\t", sum(enrich.table$positive), "/", nb.tests, " positive tests"), 3)
    }
  }

  ## Plot the different multiple testing corrections
  if (plot.adjust) {
    threshold.colors <- c("qvalue"="darkgreen", "evalue"="blue", "positive"="red")
    plot(enrich.table[,c("pvalue","evalue")], log="xy", panel.first=grid(lty="solid", col="#BBBBBB"),
         xlab="p-value", ylab="Multiple testing corrections",
         ylim=c(min(enrich.table$pvalue), max(enrich.table$evalue)),
         pch=1, col=threshold.colors["evalue"])
    abline(h=1, lwd=1)
    abline(h=thresholds, col=threshold.colors[names(thresholds)], lwd=1)  ## Mark thresholds
    abline(a=0, b=1)
    points(enrich.table[,c("pvalue","qvalue")], pch=20,
           col=threshold.colors["qvalue"])
    points(enrich.table[enrich.table["positive"]==1,c("pvalue","qvalue")], pch=3,
           col=threshold.colors["positive"], lwd=2)

    ## Count number of selected genes by different criteria
    selection.columns <- paste(sep="_", names(thresholds), thresholds)
    selection.columns <- append(selection.columns, "positive")
    selected.numbers <- apply(enrich.table[,selection.columns], 2, sum)
    legend("topleft",
           legend=paste(selected.numbers, selection.columns),
           col=threshold.colors[c(names(thresholds), "positive")],
           pch=c(1, 20, 3))
  }
  ## Select only the GO classes passing the e-value
  if (select.positives) {
    enrich.table <- enrich.table[enrich.table$positive, ]
    verbose(paste("\tGO over-representation. Returning",
                  sum(enrich.table$evalue.selection), "significant associations."), 2)
  } else {
    verbose(paste(sep="", "\tGO over-representation. Returning table with all associations (",
                  nrow(enrich.table), " among which ",
                  sum(enrich.table$positive), " significant)."), 2)
  }

  return(enrich.table)
}
