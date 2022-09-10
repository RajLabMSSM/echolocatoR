#' P-values to symbols
#'
#' Convert p-values to significance symbols.
#' @keywords internal
#' @importFrom stats symnum
pvalues_to_symbols <- function(pval_vector){
  symnum.args <- list(cutpoints = c(0,0.00001, 0.0001, 0.001, 0.01, 0.05, 1),
                      symbols = c("*****","****", "***", "**", "*", "ns"))
  symnum.args$x <- as.numeric(pval_vector)
  pvalue.signif <- do.call(stats::symnum, symnum.args) |>
    as.character()
  return(pvalue.signif)
}
