#' Extract the locus dir
#'
#' @family directory functions
#' @keywords internal
get_locus_dir <- function(subset_path){
  locus_dir <- dirname(subset_path)
  return(locus_dir)
}
