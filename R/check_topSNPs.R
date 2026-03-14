#' Check \code{topSNPs}
#'
#' Check \code{topSNPs}. Creates a data.frame from the full summary stats
#' if \code{topSNPs} is set to "auto".
#' @inheritParams finemap_loci
#' @inheritDotParams echodata::import_topSNPs
#' @returns data.frame
#'
#' @keywords internal
#' @importFrom echodata import_topSNPs
check_topSNPs <- function(topSNPs,
                          fullSS_path,
                          verbose=TRUE,
                          ...){
  topSS <- if(is.character(topSNPs) &&
              tolower(topSNPs)=="auto"){
    messager("Creating topSNPs from full summary stats file.",v=verbose)
    fullSS_path
  } else {
    topSNPs
  }
  topSNPs <- echodata::import_topSNPs(topSS = topSS,
                                      munge = TRUE,
                                      verbose = verbose,
                                      ...)
  #### Validate required columns ####
  required <- c("Locus","CHR","POS")
  ## Also accept BP as alias for POS
  if(!"POS" %in% colnames(topSNPs) && "BP" %in% colnames(topSNPs)){
      data.table::setnames(topSNPs, "BP", "POS")
  }
  missing_cols <- setdiff(required, colnames(topSNPs))
  if(length(missing_cols) > 0){
      stop("topSNPs is missing required column(s): ",
           paste(missing_cols, collapse=", "),
           "\nRequired columns: Locus, CHR, POS",
           "\nOptional columns: SNP (lead SNP RSID), Gene",
           "\nFound columns: ",
           paste(colnames(topSNPs), collapse=", "))
  }
  return(topSNPs)
}
