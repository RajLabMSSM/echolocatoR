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
  return(topSNPs)
}
