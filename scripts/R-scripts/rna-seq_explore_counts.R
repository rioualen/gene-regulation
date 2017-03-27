#' @title Descriptive statistics from a read count table.
#' @author Jacques van Helden
#' @description Compute descriptive statistics for each sample of a 
#' RNA-seq count table. The table contains one row per feature 
#' (gene, transcript) and one column per sample. 
#' @param config.file a yaml-formatted file describing the parameters 
#' of the workflow. 
#' @param count.file file containing the counts of reads per feature, 
#' tab-separated value text file.
#' @export
RnaSeqExploreCounts <- function(config.file, count.file) {
  
  ## Initialize the result object
  result <- list()
  result$config.file <- config.file
  result$count.file <- count.file
  
  ## Read the configuration file
  library(yaml)
  param <- yaml.load_file(config.file)
  
  
  ################################################################
  ## Read the count table
  verbose("Loading count table", 2)
  counts <- read.delim(count.file, row.names=1, sep="\t")

  ## Filter out the rows corresponding to non-assigned counts, 
  ## e.g. __not_aligned, __ambiguous, __too_low_qAual, __not_aligned
  not.feature <- grep(rownames(counts), pattern = "^__")
  if (length(not.feature) > 0) {
    counts <- counts[-not.feature,]
  }
  
  ## Check that some required parameters are well defined
  ## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
  epsilon <- as.logical(param$deg$epsilon)
  if (!is.numeric(epsilon)) {
    param$deg$epsilon <- 0.1
  }
  
  
  ## Default method for the selection of the final list of DEG
  selectionCriterion <- param$deg$selectionCriterion
  if (is.null(selectionCriterion)) {
    param$deg$selectionCriterion <- "edgeR"
  }
  
  ## Check that selection criterion is valid
  acceptedSelectionCriteria <- c("edgeR", "DESeq2", "union", "intersection")
  if (!selectionCriterion %in% acceptedSelectionCriteria) {
    stop("Invalid selection criterion. Supported: ", 
            paste(acceptedSelectionCriteria, collapse=", "), ".")
  }
  
  gtf.file <- param$genome$gtf_file
  
  ################################################################
  ## Load gene information from the GTF file
  ## (should be the same as used to count tags per gene)
  if (!is.null(gtf.file)) {
    ## Jacques or Lucie, TO DO : check if the method makeTxDbFromGFF allows 
    ## to get gene names and descriptions.
    ## 
    ## In principle it should be possible since this info is in the GTF file. 
    ## If not, we migh rewite a parser for GTF. 
    
    message("Loading gene annotations from GTF file: ", gtf.file)
    gtf.file.path <- file.path(dir.main, param$dir$genome, gtf.file)
    library(rtracklayer)  # for the import() function
    #   gr <- import(gtf.file)
    txdb <- makeTxDbFromGFF(file=gtf.file.path)
    #organism=organism.name, ## OPTION organism DOES NOT WORK WITH DESULFO because not supported in this library. However the GTF is correctly loaded
    #                          dataSource=gtf.source)
    #seqlevels(txdb)
    #   transcripts <- transcriptsBy(txdb, by="gene")
    #   transcript.table <- as.data.frame(transcripts@unlistData)
    #   cds <- cdsBy(txdb, by="gene")
    #   cds.table <- as.data.frame(cds@unlistData)
    #   exons <- exonsBy(txdb, by="gene")
    #   exon.table <- as.data.frame(exons@unlistData)
    all.genes <- genes(txdb)
    gene.info <- as.data.frame(all.genes)
    gene.info$name <- gene.info$gene_id
    #  gene.info$entrez.id <- NA
    gene.info$description <- "no description"
  } else {
    verbose(paste("No GTF file has been specified"))
    g <- nrow(counts)
    gene.info <- data.frame("seqnames"=rep(NA, times=g),
                            "start"=rep(NA, times=g),
                            "end"=rep(NA, times=g),
                            "width"=rep(NA, times=g),
                            "strand"=rep(NA, times=g),
                            "gene_id"=row.names(counts),
                            "name"=row.names(counts),
                            #                          "entrez.id" = rep(NA, times=g),
                            "description"=rep("no description", times=g))
    row.names(gene.info) <- row.names(counts)
  }
    
  ## Logical variable indicating whether or not the data tables should be exported 
  ## in Excel format.
  export.excel.files <- as.logical(param$deg$export_xlsx)
  if (is.na(export.excel.files)) {
    param$deg$export_xlsx <- FALSE
  }

  ## Output directory for the differential expression analysis
  dir.diffexpr <- param$dir$diffexpr
  if (is.null(dir.diffexpr))  {
    param$dir$diffexpr <- file.path("results", "diffexpr")
  }
    
  result$param <- param ## Include the parameters in the result
  result$counts <- counts
  result$gene.info <- gene.info
  return(result)
}