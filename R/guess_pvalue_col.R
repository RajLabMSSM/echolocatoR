guess_pvalue_col <- function(dat,
                             QTL_prefix=NULL){
  p_options <- c("p","p-value","p-values","pvalue","pvalues","pval")
  p_options <- c(p_options, stringr::str_to_sentence(p_options))
  pval_col <- paste0(QTL_prefix,p_options)[paste0(QTL_prefix,p_options) %in% colnames(dat)][1]
  return(pval_col)
}
