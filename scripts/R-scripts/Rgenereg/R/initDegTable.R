################################################################
## Initiate a result table with the CPMs and derived statistics
initDegTable <- function(cmps) {
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
