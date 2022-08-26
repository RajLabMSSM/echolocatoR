#' Construct the path of the locus subset
#'
#' @family directory functions
#' @keywords internal
#' @examples
#' subset_path <- echolocatoR:::construct_subset_path(dataset_type="GWAS",
#'                                                    dataset_name="Nalls2019",
#'                                                    locus="BST1")
construct_subset_path <- function(subset_path="auto",
                                  results_dir=tempdir(),
                                  dataset_type,
                                  dataset_name,
                                  locus=NULL,
                                  suffix=".tsv.gz"){
  # Specify subset file name
  if(subset_path=="auto"){
    locus_dir <- construct_locus_dir(results_dir = results_dir,
                                dataset_type = dataset_type,
                                dataset_name = dataset_name,
                                locus = locus)
    created_sub_path <- file.path(
      locus_dir, paste0(locus,"_",dataset_name,"_subset",suffix) )
    return(created_sub_path)
  } else{return(subset_path)}
}
