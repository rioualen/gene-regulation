#' @title check required libraries for DEG analysis, and install them if required
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @param required.libraries a vector contianing the names of required CRAN libraries, which will be installed with install.packages()
#' @param required.bioconductor a vector containing the required BioConductor libraries, which will be installed with biocLite
CheckRequiredLibraries <- function (required.libraries, required.bioconductor=NULL) {
  for (lib in required.libraries) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib)
      library(lib, character.only = TRUE)
    }
  }
  
  for (lib in required.bioconductor) {
    if (!require(lib, character.only = TRUE)) {
      ## try http:// if https:// URLs are not supported
      source("https://bioconductor.org/biocLite.R")
      biocLite(lib)
    }
    if (!require(lib, character.only = TRUE)) {
      stop("Missing library: ", lib, " could not be installed")
    }
  }
  
}

#' @title Load parameters and required libraries
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
LoadDEGparam <- function (yamlFile) {
  library(yaml)
  data -> yaml.load_file("vn.yamloo")
}


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
                               thresholds=thresholds,
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
  ## Beware: for the fold-change we ake the absolute value of log2FC, 
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
                                  thresholds = thresholds,
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
    if (is.na(organism.names["clusterProfiler"])) {
      message("Skipping clusterProfiler analysis because organism.names[\"clusterProfiler\"] is not specified")
    } else {
      ## Compute the enrichment of the gene list for Biological Processes of the Gene Ontology
      ## We intently set all cutoffs to a value > than the possible max, in order to get the full table.
      ## We then perform the filter only if requested.
      ego <- enrichGO(gene = geneset, 
                      organism = organism.names["clusterProfiler"],
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

#' @title Draw an equivalent of the t-test volcano plot for RNA-seq data.  
#' 
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Draw an equivalent of the t-test volcano plot for RNA-seq 
#' data. The main idea is to display simultaneously the effect size (X axis) 
#' and the significance (Y axis). The measure of the effect size is in this 
#' case the log2 fold change. (logFC). 
#' 
#' There is a close relationship between this RNA-seq version of the volcano 
#' and the volcano plots commonly used to show the result of a t-test with microarrays.
#'
#' @details
#' First version: 2015-09 
#' Last modification: 2015-09
#'
#' @param deg.table A table with the results of differential expression. Such tables can be produced by DESeq2, edgeR or other packages.
#' @param alpha=0.05 Threshold to declare features (genes) as positive or negative.
#' @param effect.size.col="log2FC" Number or name of the column containing the effect size, whose values will be displayed on the X axis. 
#' @param control.type = "p.value" Control type. Can be a p-value or an equivalent estimation of significance (FDR, e-value, FWER, ...). This value will be converted to y=-log10(x) to obtain the values displayed on the Y axis. The alpha threshold applies on the non-log-converted values.
#' @param ... 
#'
#' @examples
#'
#' @export
rnaseq.volcanoPlot <- function(deg.table,
                               alpha=0.05,
                               effect.size.col = "log2FC",
                               control.type = "p.value",
                               ...) {
  stop("NOT IMPLEMENTED YET, ASK Jacques.van-Helden@univ-amu.fr")
  volcanoPlot()
}


################################################################
## Compute sample-wise statistics
calc.stats.per.sample <- function(sample.descriptions, 
                                  count.table) {
  verbose("Computing statistics per sample", 2)
  
  stats.per.sample <- cbind(  
    sample.desc[names(count.table), ],
    data.frame(
      "zeros" = apply(all.counts == 0, 2, sum, na.rm=TRUE), ## Number of genes with 0 counts
      "detected" = apply(all.counts > 0, 2, sum, na.rm=TRUE), ## Number of genes counted at least once
      "sum" = apply(all.counts, 2, sum, na.rm=TRUE), ## Sum of all counts for the sample
      "mean" = apply(all.counts, 2, mean, na.rm=TRUE), ## Mean counts per gene
      "min" = apply(all.counts, 2, min, na.rm=TRUE), ## Min counts per gene 
      "perc05" = apply(all.counts, 2, quantile, probs=0.05, na.rm=TRUE), ## 5th percentile
      "perc25" = apply(all.counts, 2, quantile, probs=0.25, na.rm=TRUE), ## 25th percentile
      "median" = apply(all.counts, 2, median, na.rm=TRUE), ## median (percentile 50)
      "perc75" = apply(all.counts, 2, quantile, probs=0.75, na.rm=TRUE), ## percentile 75
      "perc95" = apply(all.counts, 2, quantile, probs=0.95, na.rm=TRUE), ## percentile 95
      "max" = apply(all.counts, 2, max, na.rm=TRUE) ## Max counts per gene
    )
  )
  stats.per.sample$max.sum.ratio <- stats.per.sample$max / stats.per.sample$sum
  stats.per.sample$median.mean.ratio <- stats.per.sample$media / stats.per.sample$mean
  stats.per.sample$Mreads <- round(stats.per.sample$sum/1e6, digits = 1)

  ## Count number and the fraction of samples with counts below the mean. 
  ## This shows the impact of very large counts: in some samples, 
  ## 85% of the samples have a value below the mean (i.e. the mean is at the percentile 85 !)
  stats.per.sample$below.mean <- apply(t(all.counts) < stats.per.sample$mean, 1, sum, na.rm=TRUE)
  stats.per.sample$fract.below.mean <- stats.per.sample$below.mean/nrow(all.counts)
  # View(stats.per.sample)
  
  return(stats.per.sample)
}

################################################################
## Draw a barplot with the number of reads per sample
libsize.barplot <- function(stats.per.sample, 
                            main="Read library sizes (libsum per sample)",
                            plot.file=NULL, 
                            ...) {
  
  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- max(nchar(sample.desc$label))/3+5
  boxplot.height <- length(sample.ids)/3+2
  
  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }
  
  par(mar=c(5,boxplot.lmargin,4,1)) ## adapt axes
  bplt <- barplot(stats.per.sample$Mreads, names.arg = stats.per.sample$label, 
                  main=main,
                  horiz = TRUE, las=1,
                  xlab="libsum (Million reads per sample)",
                  col=stats.per.sample$color, ...)
  grid(col="white", lty="solid",ny = 0)
  text(x=pmax(stats.per.sample$Mreads, 3), labels=stats.per.sample$Mreads, y=bplt,pos=2, font=2)
  if (!is.null(plot.file)) {
    silence <- dev.off()
  }
}

################################################################
## Draw boxplots with read counts per genes for each sample
count.boxplot <- function(count.table, 
                          sample.desc,
                          sample.label.col=1,
                          xlab="Raw counts",
                          main="Box plots per sample: raw counts",
                          plot.file=NULL) {
  
  ## Adapt boxplot size to the number of samples and label sizes
  boxplot.lmargin <- max(nchar(sample.desc$label))/3+5
  boxplot.height <- length(sample.ids)/3+2
  
  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }
  
  par(mar=c(5,boxplot.lmargin,4,1)) ## adapt axes
  boxplot(count.table, horizontal=TRUE, col=sample.desc$color,
          xlab=xlab, names=sample.desc[, sample.label.col],
          main=main, las=1)
  grid(col="grey", lty="solid",ny = 0)
  if (!is.null(plot.file)) {
    silence <- dev.off()
  }
}

################################################################
## Draw a heatmap with the inter-sample correlation matrix.
count.correl.heatmap <- function(count.table, 
                                 main="Correlation between raw counts",
                                 plot.file=NULL,
                                 log.transform=FALSE, # Perform a log transformation of the values before plotting
                                 epsilon=0.1, # Add an epsilon to zero values before log transformation, in order to -Inf values
                                 ...
                                 ) {
  
  
  ## Adapt boxplot size to the number of samples and label sizes
  margin <- max(nchar(names(count.table)))/3+5

  if (log.transform) {
    range(count.table)
    count.table[count.table==0] <- epsilon
    count.table <- log10(count.table)
  }
  count.cor <- as.matrix(cor(count.table))
  
  ## Define a color palette for heatmaps. I like this Red-Blue palette because 
  ## - it suggests a subjective feeling of warm (high correlation)/cold (low correlation)
  ## - it can be seen by people suffering from red–green color blindness.
  cols.heatmap <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))
  
  ## Use a grayscale color  
#
  cols.heatmap <- gray.colors(100, start = 1, end = 0, gamma = 3, alpha = NULL)

  ## Sample-wise library sizes
  if (!is.null(plot.file)) {
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=8, height=boxplot.height)
  }
  
  hm <- heatmap.2(count.cor,  scale="none", trace="none", 
                  #breaks=c(-1, seq(0,1,length.out = 100)),
                  main=main, margins=c(margin,margin),
                  col=cols.heatmap,
                  cellnote = signif(digits=2, count.cor),
                  ...
                  )
  

  if (!is.null(plot.file)) {
    silence <- dev.off()
  }

  return(count.cor)  
}

## Generate a set of plots displaying some sample-wise statistics
sample.description.plots <- function (sample.desc,
                                      stats.per.sample,
                                      dir.figures,
                                      exploratory.plots=FALSE) {
  
  plot.files <- c() ## To retun: list of files with plots
  
  par.ori <- par() ## Save original plot parameters
  
  
  ## Library size barplots
  plot.files["Mreads_barplot"] <- file.path(dir.figures, "sample_libsum_barplot.pdf")
  libsize.barplot(stats.per.sample, plot.files["Mreads_barplot"])
  
  ################################################################
  ## Boxplots of raw counts and CPMs, in linear + log axes. 
  ## These plots give a pretty good intuition of the raw data per sample: 
  ## library sizes, outliers, dispersion of gene counts.
  
  ## Boxplot of raw counts
  plot.files["sample_boxplot_counts"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_counts.pdf"))
  count.boxplot(all.counts, stats.per.sample, 
                xlab="Raw counts", main="Box plots per sample: raw counts",
                plot.file=plot.files["sample_boxplot_counts"])

  ## Boxplot of log10-transformed counts
  plot.files["sample_boxplot_counts_log10"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_counts_log10.pdf"))
  count.boxplot(all.counts.log10, stats.per.sample, 
                xlab="log10(counts)", main="Box plots per sample: log10(counts)",
                plot.file=plot.files["sample_boxplot_counts_log10"])

  ## Boxplot of CPMs
  plot.files["sample_boxplot_CPM"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_CPM.pdf"))
  count.boxplot(cpms, stats.per.sample, 
                xlab="CPM", main="Box plots per sample: counts per million reads (CPM)",
                plot.file=plot.files["sample_boxplot_CPM"])
  
  ## Boxplot of log10-transformed CPMs
  plot.files["sample_boxplot_CPM_log10"] <- file.path(dir.figures, paste(sep = "", "sample_boxplots_CPM_log10.pdf"))
  count.boxplot(cpms.log10, stats.per.sample, 
                xlab="log10(CPM)", main="Box plots per sample: counts per million reads (CPM)",
                plot.file=plot.files["sample_boxplot_CPM"])

  par <- par.ori ## Restore original plot parameters
  par(mar=c(4.1,5.1,4.1,1.1))

  ## Draw sample correlation heatmaps for the raw read counts
  plot.files["sample_correl_heatmap_counts"] <- paste(sep="", prefix["general.file"],"_sample_correl_heatmap_counts.pdf")
  count.correl.heatmap(all.counts, plot.file=plot.files["sample_correl_heatmap_counts"])
#   hm <- heatmap.2(,  scale="none", trace="none", 
#                   main="Correlation between raw counts", margins=c(8,8),
#                   col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))

  ## Draw sample correlation heatmaps for CPM. Actually it gives exactly the 
  ## same result as correlation between raw counts, since the correlation has a 
  ## standardizing effect. 
  # pdf(file=paste(sep="", prefix["general.file"],"_sample_correl_heatmap_cpms.pdf"))
  # hm <- heatmap.2(as.matrix(cor(cpms)),  scale="none", trace="none", 
  #                 main="Correlation between CPM",
  #                 col=cols.heatmap) #, breaks=seq(-1,1,2/length(cols.heatmap)))
  # quiet <- dev.off()
  
  ## Plot the first versus second components of samples
  cpms.pc <- prcomp(t(cpms))
  plot.file <- paste(sep="", prefix["general.file"],"_CPM_PC1-PC2.pdf")
  plot.files["CPM_PC1-PC2"] <- plot.file
  message("Generating plot", plot.file)
  pdf(file=plot.file)
  plot(cpms.pc$x[,1:2], panel.first=grid(), type="n", main="First components from PCA-transformed CPMs")
  text(cpms.pc$x[,1:2], labels = sample.conditions, col=sample.desc$color)
  quiet <- dev.off()
  
  
  
  ## Exploratory plots, should not be done for all projects.
  if (exploratory.plots) {
    verbose("Drawing generic plots from the whole count table", 1)
    
    ## Plot the impact of the normalization factor (library sum , median or percentile 75)
    plot.file <- file.path(dir.DEG, paste(sep = "", "CPM_libsum_vs_median_vs_perc75.png"))
    plot.files["CPM_libsum_vs_median_vs_perc75"] <- plot.file
    message("Generating plot", plot.file)
    png(file= plot.file, width=1000, height=1000)
    cols.counts <- as.data.frame(matrix(sample.desc$color, nrow=nrow(all.counts), ncol=ncol(all.counts), byrow = TRUE))
    colnames(cols.counts) <- names(all.counts)
    rownames(cols.counts) <- rownames(all.counts)
    plot(data.frame("libsum" = as.vector(as.matrix(cpms.libsum)),
                    "median" = as.vector(as.matrix(cpms.median)),
                    "perc75" = as.vector(as.matrix(cpms.perc75))),
         col=as.vector(as.matrix(cols.counts)))
    quiet <- dev.off()
    
    ## Plot some sample-wise statistics
    plot.file <- file.path(dir.DEG, paste(sep = "", "sample_statistics_plots.pdf"))
    plot.files["sample_statistics_plots"] <- plot.file
    message("Generating plot", plot.file)
    pdf(file=plot.file, width=10, height=10)
    par(mar=c(5,5,1,1)) ## adpt axes
    par(mfrow=c(2,2))
    ## Median versus mean
    plot(stats.per.sample[,c("mean", "median")], 
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)
    
    ## First versus third quartile
    plot(stats.per.sample[,c("perc25", "perc75")], 
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)
    
    ## Sum versus third quartile. 
    plot(stats.per.sample[,c("sum", "perc75")], 
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)
    
    ## Mean versus third quartile. 
    plot(stats.per.sample[,c("mean", "perc75")], 
         panel.first=c(grid(lty="solid", col="#DDDDDD"), abline(a=0,b=1)),
         las=1, col=sample.desc$color)
    par(mfrow=c(1,1))
    quiet <- dev.off()
  }
  
  return(plot.files)
}

################################################################
## Initiate a result table with the CPMs and derived statistics
init.deg.table <- function(cmps) {
  verbose("\t\tInitializing result table for one differential analysis (two-sample comparison).", 2)
  all.gene.ids <- row.names(cpms)
  result.table <- data.frame("gene_id" = all.gene.ids,
                             "name"=gene.info[all.gene.ids,"name"])
  row.names(result.table) <- all.gene.ids
  result.table$entrez.id <- gene.info[all.gene.ids,"entrez.id"]
  result.table$description <- gene.info[all.gene.ids,"description"]
  
  result.table <- cbind(result.table, counts=current.counts) ## Include the original counts in the big result table
  result.table <- cbind(result.table, cpm=current.cpms) ## Include CPMs in the big result table
  
  ## Tag genes detected in less than min.rep samples, which is defined as 
  ## the minimal number of replicates per condition.
  min.rep <- min(length(samples1), length(samples2))
  result.table$undetected <- rowSums(current.counts > 1) < min.rep
  # table(result.table$undetected)
  # dim(current.counts)
  
  result.table$cpm.mean <- apply(cpms,1, mean)
  result.table$cpm1.mean <- apply(as.data.frame(cpms[,samples1]),1, mean)
  result.table$cpm2.mean <- apply(as.data.frame(cpms[,samples2]),1, mean)
  result.table$A = log2(result.table$cpm1.mean*result.table$cpm2.mean)/2
  result.table$M = log2(result.table$cpm1.mean/result.table$cpm2.mean)
  result.table$cpm.median <- apply(cpms,1, median)
  result.table$cpm1.median <- apply(as.data.frame(cpms[,samples1]),1, median)
  result.table$cpm2.median <- apply(as.data.frame(cpms[,samples2]),1, median)
  result.table$cpm.min <-  apply(cpms,1, min)
  result.table$cpm1.min <- apply(as.data.frame(cpms[,samples1]),1, min)
  result.table$cpm2.min <- apply(as.data.frame(cpms[,samples2]),1, min)
  result.table$cpm.max <-  apply(cpms,1, max)
  result.table$cpm1.max <- apply(as.data.frame(cpms[,samples1]),1, max)
  result.table$cpm2.max <- apply(as.data.frame(cpms[,samples2]),1, max)
  result.table$cpm.sd <-  apply(cpms,1, sd)
  result.table$cpm1.sd <- apply(as.data.frame(cpms[,samples1]),1, sd)
  result.table$cpm2.sd <- apply(as.data.frame(cpms[,samples2]),1, sd)
  result.table$cpm.var <-  apply(cpms,1, var)
  result.table$cpm1.var <- apply(as.data.frame(cpms[,samples1]),1, sd)
  result.table$cpm2.var <- apply(as.data.frame(cpms[,samples2]),1, sd)
  return(result.table)
}

################################################################
## DESeq2 analysis
## 
## Detect differentially expressed genes (DEG) using the package DESeq2, 
## and add a few custom columns (e-value, ...) to the result table. 
deseq2.analysis <- function(dir.figures=NULL) {
  verbose("\t\tDESeq2 analysis", 2)
  
  ## Path prefix to save DESeq2 result files
  prefix["DESeq2_file"] <- paste(sep="", prefix["comparison_file"], "_", suffix.DESeq2)
  prefix["DESeq2_figure"] <- paste(sep="", prefix["comparison_figure"], "_", suffix.DESeq2)
  
  ## Create a DESeqDataSet object from the count table + conditions
  condition <- as.factor(as.vector(current.sample.conditions))
  deseq2.dds <- DESeqDataSetFromMatrix(
    countData = current.counts, 
    colData = DataFrame(condition),
    ~ condition)
  
  
  ## Indicate that second condition is the reference condition. 
  ## If not done, the conditions are considered by alphabetical order, 
  ## which may be misleading to interpret the log2 fold changes. 
  deseq2.dds$condition <- relevel(deseq2.dds$condition, ref=cond2) 
  
  ## Run the differential analysis
  deseq2.dds <- DESeq(deseq2.dds)      ## Differential analysis with negbin distrib
  deseq2.res <- results(deseq2.dds, independentFiltering=FALSE, pAdjustMethod = "BH")  ## Collect the result table
  
  deseq2.result.table <- data.frame(
    "gene.id" = row.names(deseq2.res),
    "mean" = deseq2.res$baseMean,
    "log2FC" = deseq2.res$log2FoldChange,
    "pvalue" = deseq2.res$pvalue,
    "padj" = deseq2.res$padj)
  deseq2.result.table <- complete.deg.table(
    deseq2.result.table, 
    paste(sep="_", "DESeq2", prefix["comparison"]),
    sort.column = "padj",
    thresholds=thresholds,
    round.digits = 3,
    dir.figures=dir.figures)
  return(deseq2.result.table)
}  

edger.analysis <- function(dir.figures=NULL) {
  verbose("\t\tedgeR analysis", 2)
  
  prefix["edgeR_file"] <- paste(sep="", prefix["comparison_file"], "_", suffix.edgeR)
  prefix["edgeR_figure"] <- paste(sep="", prefix["comparison_figure"], "_", suffix.edgeR)
  
  ## Convert the count table in a DGEList structure and compute its parameters.
  d <- DGEList(counts=current.counts, group=sample.conditions[names(current.counts)])
  d$samples$group <- relevel(d$samples$group, ref=cond2) ## Ensure that condition 2 is considered as the reference
  d <- calcNormFactors(d, method="RLE")                 ## Compute normalizing factors
  d <- estimateCommonDisp(d, verbose=FALSE)             ## Estimate common dispersion
  d <- estimateTagwiseDisp(d, verbose=FALSE)            ## Estimate tagwise dispersion
  
  ################################################################
  ## Detect differentially expressed genes by applying the exact 
  ## negative binomial test from edgeR package.
  edger.de <- exactTest(d, pair=c(cond2, cond1))      ## Run the exact negative binomial test
  
  ## Sort genes by increasing p-values, i.e. by decreasing significance
  edger.tt <- topTags(edger.de, n=nrow(d), sort.by = "PValue")
  
  ## Complete the analysis of edgeR result table
  edger.result.table <- data.frame("gene.id" = row.names(edger.tt$table),
                                   "mean"=edger.tt$table$logCPM,
                                   "log2FC"=edger.tt$table$logFC,
                                   "pvalue"=edger.tt$table$PValue,
                                   "padj"=edger.tt$table$FDR)
  edger.result.table <- complete.deg.table(
    deg.table = edger.result.table, 
    table.name = paste(sep="_","edgeR",prefix["comparison"]),
    sort.column = "padj",
    thresholds=thresholds,
    round.digits = 3,
    dir.figures=dir.figures)
  # dim(edger.result.table)
  # dim(edger.result.table)
  # names(edger.result.table)
  # View(edger.result.table)
  return (edger.result.table)
}
