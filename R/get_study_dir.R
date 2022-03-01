#' Extract the study dir
#'
#' @family directory functions
#' @keywords internal
get_study_dir <- function(locus_dir){
  study_dir <- dirname(locus_dir)
  return(study_dir)
}
