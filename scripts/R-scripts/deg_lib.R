
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
  ## Select genes passing all thresholds
  deg.table[,"DEG"] <- 
    1*(apply(deg.table[,selection.columns],1,sum) == length(thresholds))
  
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

#' @title Run functional enrichment for a given gene set
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Run functional enrichment for a given gene set. 
#' This function is a compilation of pieces of code from the vignettes 
#' of different functional enrichment analysis package.
#'
#' @details
#' First version: 2015-08 
#' Last modification: 2015-08. 
#'
#' @param geneset   Entrez IDs of the gene set of interest.
#' @param allgenes  EntrezIDs for all genes of the analyzed genome (the "universe").
#' @param db an annotation database
#' @param thresholds["evalue"]=0.1 Upper threshold on the e-value for GOstats. 
#' Note that we perform the selection with a p-value threshold of 1, in order to obtain the 
#' stats for all the classes. The e-value is calculated a posteriori and indicated in the table.
#' @param select.positives=TRUE If TRUE, the result table only contains the positive tests, 
#' i.e. those passing all the thresholds. Otherwise, all rows are returned, with a column "positive" 
#' indicating whether or not they did or not pass all thresholds. Additional columns indicate the 
#' result of the individual thresholds.
#' @param organism organism name in various tools. E.g: 
#'
#' @examples
#'
#' @export
functional.enrichment <- function(geneset,
                                  allgenes,
                                  db,
                                  organism.names,
                                  ontology="BP",
                                  thresholds = c("evalue"=1, "qvalue"=0.05),
                                  select.positives=TRUE,
                                  run.GOstats = TRUE,
                                  run.clusterProfiler = FALSE,
                                  plot.adjust = TRUE) {
  #library("ALL")
  verbose(paste(sep="", "Go over-representation analysis. ",
                length(geneset), " input genes (among ", length(allgenes), ")"), 2)
  
  ## Prepare the result
  result <- list()
  result$ngenes.test <- length(geneset)
  result$ngenes.all <- length(allgenes)
  result$db = db
  result$ontology = ontology
  result$thresholds <- thresholds
  result$select.positives <- select.positives
  
  ## Create a data frame with results per gene (one row per gene, one column per attribute)
  result.per.gene <- data.frame("entrez.id" = allgenes)
  row.names(result.per.gene) <- allgenes
  result.per.gene$selected <- 0
  result.per.gene[geneset,"selected"] <- 1
  
  ## Load annotation library
  library("GO.db")
  library("annotate")
  library(db, character.only = TRUE)
  envir <- c(
    db=db,
    prefix = sub(db, pattern = ".db", replacement = ""))
  for (suffix in c("GO", "ENTREZID")) {
    envir[suffix] <- paste(sep="", envir["prefix"],suffix)
  }
  
  ################################################################
  ## Run over-representation analysis with GOstats
  if (run.GOstats) {
    library("GOstats")
    
    ## Get go annotations per gene
    go.annot.geneset <- mget(geneset, envir=get(envir["GO"]))
    go.annot.allgenes <- mget(allgenes, envir=get(envir["GO"]))
    result.per.gene$nb.go.annot <- as.vector(unlist(lapply(go.annot.allgenes, length)  ))
    
    ## Compute number of genes with at least one annotation in GO. 
    ## Note that this is for all the GO ontologies. For the hypergeometric test 
    ## below we need to take into account the ontology-specific annotations.
    nb.allgenes.with.go <- sum(result.per.gene$nb.go.annot >0) ## Total number of genes with GO annotations
    nb.geneset.with.go <- sum(result.per.gene[geneset, "nb.go.annot"] >0) ## Number of genes with GO annotations in the test set
    # head(result.per.gene)
    # table(result.per.gene$nb.go.annot)
    
    ## Define parameters for GOstats analysis
    gostats.params <- new("GOHyperGParams",
                          geneIds=geneset,
                          universeGeneIds=allgenes,
                          annotation=db,
                          ontology=ontology,
                          pvalueCutoff=1,
                          conditional=FALSE,
                          testDirection="over")

    ## Run over-representation hypergeometric test with GOstats
    over.result <- hyperGTest(gostats.params)
    # plot(goDag(over.result))
    
    result$geneIds <-  length(geneIds(over.result)) # Number of mapped genes
    result$geneIdUniverse <- length(geneIdUniverse(over.result))
    
    ## Collect the statistics for all the tested GO classes.
    ## For this we set the p-value threshold to 1.1, so we select
    ## all classes, even non-significant!
    go.enrich.table <- summary(over.result, pvalue=2)
    go.enrich.table$expectedCounts = expectedCounts(over.result)
    
    go.enrich.table <- complete.enrich.table(go.enrich.table,
                                             pvalue.column="Pvalue",
                                             thresholds=thresholds, 
                                             select.positives=select.positives,
                                             plot.adjust=plot.adjust)    
    
    ## Add  GO table to the result
    result$go.bp.table <- go.enrich.table
    
    
    ## For the sake of understanding, we can re-compute the p-value 
    ## with the hypergeometric distribution. We obtain the same result as GOstats.
    N <- universeMappedCount(over.result) ## Universe, i.e. nb of genes with at least one annotation in the selected ontology
    m <- go.enrich.table$Size ## genes marked as belonging to the considered GO class
    n <- N -m ## Number of "non-marked" genes, i.e. not belonging to the considered GO class
    k <- length(over.result@geneIds) ## Number of test genes with at least one annotation in the selected ontology
    Pvalue.check <- 
      phyper(q = go.enrich.table$Count -1, m=m, n=n, k = k, lower=FALSE)
    #   plot(go.enrich.table[,pvalue.column], Pvalue.check, log="xy")
    ExpCount.check <- m * k/N
    # plot(go.enrich.table$ExpCount, ExpCount.check)
    #   abline(a=0, b=1)
    
    # View(go.enrich.table)
    # names(go.enrich.table)
    # dim(go.enrich.table)
    
  }    ## End of GOstats analysis
  
  
  ## Run GO enrichment analysis with clusterProfiler
  if (run.clusterProfiler) {
    library("clusterProfiler")
    ## Check the organism
    if (is.na(organism["clusterProfiler"])) {
      message("Skipping clusterProfiler analysis because organisms is not specified")
    } else {
      ## Compute the enrichment of the gene list for Biological Processes of the Gene Ontology
      ## We intently set all cutoffs to a value > than the possible max, in order to get the full table.
      ## We then perform the filter only if requested.
      ego <- enrichGO(gene = geneset, organism = organism["clusterProfiler"],
                      pvalueCutoff = 2,
                      qvalueCutoff = 2,
                      ont = ontology, readable = T)  
      ego.table <- data.frame(ego@result)

      
      ego.table <- complete.enrich.table(
        ego.table,
        pvalue.column="pvalue",
        thresholds=thresholds, 
        select.positives=select.positives,
        plot=TRUE)
      
      ## Add  GO table to the result
      result$clusterProfiler.bp.table <- ego.table
    }
  }
  
  ## Att the results per gene to the result object
  result$result.per.gene <- result.per.gene
  
  # plot(hist(go.enrich.table[,pvalue.column], breaks=20))
  #dim(go.enrich.table)
  # View(go.enrich.table)
  # print(over.result.cond)
  return(result)
}

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
complete.enrich.table <- function(enrich.table, 
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
