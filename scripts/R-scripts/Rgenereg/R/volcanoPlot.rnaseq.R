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
volcanoPlot.rnaseq <- function(deg.table,
                               alpha=0.05,
                               effect.size.col = "log2FC",
                               control.type = "p.value",
                               ...) {
  stop("NOT IMPLEMENTED YET, ASK Jacques.van-Helden@univ-amu.fr")
  volcanoPlot()
}
