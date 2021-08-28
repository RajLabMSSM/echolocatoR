#' @family general
#' @keywords internal
get_header <- function(large_file,
                       colnames_only=TRUE,
                       n=2,
                       nThread=1){
  ### Reading in this way is more robust and able to handle bgz format.
  header <- data.table::fread(text=readLines(con = large_file, n = n),nThread = nThread)
  if(colnames_only) header <- colnames(header)
  return(header)
}



