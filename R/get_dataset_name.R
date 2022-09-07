#' @family directory functions
#' @keywords internal
#' @importFrom utils tail
get_dataset_name <- function(file_path){
  dataset_name <- utils::tail(strsplit(dirname(file_path), "/")[[1]],n = 1)
  return(dataset_name)
}
