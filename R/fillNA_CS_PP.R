#' Fill NA in PP and CS columns
#'
#' @family general
#' @keywords internal
#' @examples
#' data("BST1");
#' finemap_dat <- BST1
#' # finemap_dat <- data.table::fread("~/Desktop/results/GWAS/Kunkle_2019.microgliaQTL/ABCA7/Multi-finemap/ABCA7_Kunkle_2019.microgliaQTL_Multi-finemap.tsv.gz")
#' finemap_dat <- fillNA_CS_PP(finemap_dat=finemap_dat)
fillNA_CS_PP <- function(finemap_dat,
                         fillNA_CS=0,
                         fillNA_PP=0){
  CS_cols <- grep(".CS$",colnames(finemap_dat), value = T)
  PP_cols <-  grep(".PP$",colnames(finemap_dat), value = T)
  finemap_dat <- data.frame(finemap_dat)
  if(!is.null(fillNA_CS)){
    printer("+ Filling NAs in CS cols with",fillNA_CS)
    finemap_dat[,CS_cols][is.na(finemap_dat[,CS_cols])] <- fillNA_CS
  }
  if(!is.null(fillNA_PP)){
    printer("+ Filling NAs in PP cols with",fillNA_PP)
    finemap_dat[,PP_cols][is.na(finemap_dat[,PP_cols])] <- fillNA_PP
  }
  return(finemap_dat)
}



