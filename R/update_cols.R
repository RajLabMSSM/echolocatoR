#' Update CS cols
#'
#' Convert old col format: ".Credible_Set" => ".CS"
#' @examples
#' data("BST1");
#' finemap_DT <- update_cols(finemap_dat=BST1)
update_cols <- function(finemap_dat){
  colnames(finemap_dat) <- gsub("*.Credible_Set$",".CS",colnames(finemap_dat))
  colnames(finemap_dat) <- gsub("*.Probability$",".PP",colnames(finemap_dat))
  return(finemap_dat)
}



