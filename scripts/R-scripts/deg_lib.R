
#' @title Display messages at a given verbosity level
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Display messages depending on user-defined verbosity level
#'
#' @details
#' First version: 2015-03. 
#' Last modification: 2015-03. 
#'
#' @param verbosity   Level of verbosity above which the message should be printed.
#' @param print.date=TRUE   Print date and time
#'
#' @examples
#'
#' verbosity <- 1 ## Define level of verbosity
#'
#' ## This message will be printed because the level is <= verbosity
#' verbose("This is printed", 1)
#'
#' ## This message will not be printed because the verbosity is inferior to the specified level
#' verbose("This is not printed", 2)
#'
#' @export
verbose <- function(message.content,
                    level=1,
                    print.date=TRUE) {
  if (!exists("verbosity")) {
    verbosity <- 1
  }
  if (verbosity >= level) {
    if (print.date) {
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t", message.content)
    } else {
      message(message.content)
    }
  }
}

#' @title Analyse a table of differentially expressed genes.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Analyse a table of differentially expressed genes to 
#' produce plots and additional statistics. 
#'
#' @details
#' First version: 2015-08 
#' Last modification: 2015-08. 
#'
#' @param deg.table   Table with the results of differential expression analysis
#' obtained from RNA-seq data. 
#' Such tables can be produced by a variety of packages, e.g. DESeq2, edgeR, etc. 
#' However, the column names should be redefined since each of this packages uses 
#' different names for similar metrics. 
#' The required columns are c("gene.id", mean", "log2FC", "pvalue", "padj").
#' Additional columns can be provided but will be ignored for the analysis.
#' 
#' @param table.name  name for the data table, which will be displayed in the plots. 
#' 
#' @param sort.column="none" Column to sort the result. 
#' Supported: c("none", "mean", "log2FC", "pvalue", "padj").
#' 
#' @param thresholds=c("pvalue" = 0.05, "padj"=0.05, "evalue"=1, "FC"=1.5) 
#' Thresholds on some specific scores, used for display purpose, and to 
#' add some columns to the result table, indicating if the gene passes the 
#' threshold or not.
#' 
#' @param round.digits=3 Significant digits to round the values of the 
#' output table. Set to NA to avoid rounding. 
#'
#' @param dir.figures=NULL if not NULL, figures will be saved in the specified directory. 
#' 
#' @examples
#'
#'  deseq2.result.table <- data.frame(
#'     "gene.id" = row.names(deseq2.res),
#'     "mean" = deseq2.res$baseMean, 
#'     "log2FC" = deseq2.res$log2FoldChange,
#'     "pvalue" = deseq2.res$pvalue,
#'     "padj" = deseq2.res$padj)
#' deseq2.result.table <- complete.deg.table(deseq2.result.table, 
#'     table.name="DESeq2", sort.column="padj")
#'
#' @export
complete.deg.table <- function(deg.table,
                               table.name,
                               sort.column = "none",
                               thresholds=c("padj"=0.05, "evalue"=1, "FC"=1.5),
                               round.digits = 3,
                               dir.figures=NULL) {
  
  # names(deg.table)
  verbose(paste("Analysing", table.name, "result table with", nrow(deg.table), "rows"), 1)

  col.descriptions <- vector() ## Initialize vector with column descriptions
  
  # Check that the input table contains the required columns
  required.columns <- c("gene.id", "mean", "log2FC", "pvalue", "padj")
  for (column in required.columns) {
    if (!column %in% names(deg.table)) {
      stop(paste("DEG table should contain",column,"column"))
    }
  }
  
  ## Make sure that row names correspond to gene IDs
  row.names(deg.table) <- deg.table$gene.id
  
  ## Sort DEG table if required
  sort.decreasing <- c("mean"=TRUE, 
                       "padj"=FALSE, 
                       "pvalue"=FALSE, 
                       "log2FC"=TRUE)
  if (sort.column != "none") {
    verbose(paste("\tSorting DEG table by", sort.column), 2)    
    deg.table <- deg.table[order(deg.table[,sort.column], 
                                 decreasing = sort.decreasing[sort.column]),]  
  }
  # head(deg.table)
  # tail(deg.table)
  # summary(deg.table)
  
  ## Compute fold-change from log2 fold change
  ## Beware: for the fold-chante we ake the absolute value of log2FC, 
  ## to have a fold change irrespective of the up- or -down sense
  deg.table$FC <- 2^abs(deg.table$log2FC) 
  
  col.descriptions["FC"] <- "Sense-insensitive fold change (always >= 1)"
  as.data.frame(col.descriptions)
  
  ## Compute E-value
  deg.table$evalue <- deg.table$pvalue * nrow(deg.table)
  col.descriptions["evalue"] <- "Expected nb of FP (=pvalue * number of tests)"
  
  ## Compute the rank of p-value sorted genes
  deg.table$padj.rank <- rank(deg.table$padj)
  col.descriptions["pval.rank"] <- "Rank of genes sorted by increasing adjusted p-values"
  
  ## Indicate the sense of the regulation (up-down)
  deg.table$sign <- sign(deg.table$log2FC)
  col.descriptions["sign"] <- "Sign of the regulation (1= up, -1 = down)"
  
  ## Label the genes passing the FDR, E-value and fold-change thresholds
  threshold.type <- c("pvalue"="upper", "padj"="upper", "evalue"="upper", "FC"="lower")
  selection.columns <- paste(sep="", names(thresholds), "_", thresholds)
  names(selection.columns) <- names(thresholds)
  for (s in names(thresholds)) {
    if (threshold.type[s] == "upper") {
      selected <- deg.table[, s] < thresholds[s]
    } else {
      selected <- deg.table[, s] > thresholds[s]
    }
    deg.table[, selection.columns[s]] <- selected*1
    col.descriptions[selection.columns[s]] <- paste("Passing", threshold.type[s], "threshold on", s)
  }
  
  # print(data.frame(col.descriptions))
  
  # table(deg.table[,selection.columns])
  
  ## Round columns to a reasonable number of significant digits
  if (!is.na(round.digits)) {
    verbose(paste("Rounding values to", round.digits, "digits"), 2)
    for (col in setdiff(names(deg.table), selection.columns)) {
      if (is.numeric(deg.table[,col])) {
        verbose(paste("Rounding", col), 2)
        deg.table[,col] <- signif(digits=3, deg.table[,col])
      }    
    }
  }
  
  ################################################################
  ## Export figures
  if (!is.null(dir.figures)) {
    verbose(paste("\tSaving figures in directory", dir.figures), 2)
    dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)
    
    ## Draw Venn diagram with number of genes declared significant 
    ## according to the selection criteria (threshold fields).
    selection.venn.counts <- vennCounts(deg.table[,selection.columns])
    pdf(file=file.path(dir.figures, paste(sep = "", table.name, "selection_Venn.pdf")))
    vennDiagram(selection.venn.counts, cex=1, main=paste(table.name, "selected genes"))  
    silence <- dev.off()
  
    ## Histogram of the nominal p-values
    pdf(file=file.path(dir.figures, paste(sep = "", table.name, "_pval_hist.pdf")), width=7, height=5)
    hist(deg.table$pvalue, breaks=seq(from=0, to=1, by=0.05),
         xlab="Nominal p-value",
         ylab="Number of genes",
         main=paste(table.name, "pvalue distribution"),
         col="#BBBBBB")
    silence <- dev.off()
  }
  return(deg.table)
}
