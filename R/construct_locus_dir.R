#' Make locus-specific results folder
#'
#' @family directory functions
#' @keywords internal
#' @examples
#' locus_dir <- echolocatoR:::construct_locus_dir(dataset_type="GWAS",
#'                                                dataset_name="Nalls2019",
#'                                                locus="BST1")
construct_locus_dir <- function(results_dir=tempdir(),
                                dataset_type,
                                dataset_name,
                                locus){
  locus_dir <- file.path(results_dir, dataset_type, dataset_name, locus)
  dir.create(locus_dir, showWarnings = FALSE, recursive = TRUE)
  return(locus_dir)
}
