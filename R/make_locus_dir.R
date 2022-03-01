#' Make locus-specific results folder
#'
#' @family directory functions
#' @keywords internal
#' @examples
#' locus_dir <- echolocatoR:::make_locus_dir(results_dir = file.path(tempdir(),
#'                                                                   "results"),
#'                                           dataset_type="GWAS",
#'                                           dataset_name="Nalls23andMe_2019",
#'                                           locus="BST1")
make_locus_dir <- function(results_dir,
                           dataset_type="dataset_type",
                           dataset_name="dataset_name",
                           locus){
  locus_dir <- file.path(results_dir, dataset_type, dataset_name, locus)
  dir.create(locus_dir, showWarnings = FALSE, recursive = TRUE)
  return(locus_dir)
}
