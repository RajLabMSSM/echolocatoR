#' @family directory functions
#' @keywords internal
get_dataset_name <- function(file_path){
  dataset_name <- tail(strsplit(dirname(file_path), "/")[[1]],n = 1)
  return(dataset_name)
}
