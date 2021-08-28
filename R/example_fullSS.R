#' Store \code{\link{fullSS_dat}}
#'
#' @inheritParams finemap_pipeline
#' @family general
#' @keywords internal
#' @examples
#' \dontrun{
#' data("fullSS_dat")
#' fullSS_path <- example_fullSS()
#' }
example_fullSS <- function(fullSS_path=file.path(tempdir(),"Nalls23andMe_2019.fullSS_subset.tsv"),
                           munged=TRUE,
                           nThread=1){
  dat <- if(munged) echolocatoR::fullSS_munged else echolocatoR::fullSS_dat
  printer("Writing file to ==>",fullSS_path)
  data.table::fwrite(x = dat,
                     file = fullSS_path,
                     nThread = nThread,
                     sep="\t")
  return(fullSS_path)
}



