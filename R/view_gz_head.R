#' @family general
#' @keywords internal
view_gz_head <- function(gz_path, nrow=10){
  system(paste("zcat",gz_path,"| head",nrows))
}



