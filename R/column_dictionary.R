#' Map column names to positions.
#'
#' Useful in situations where you need to specify columns by index instead of name (e.g. awk queries).
#'
#' @param file_path Path to full summary stats file
#' (or any really file you want to make a column dictionary for).
#' @return Named list of column positions.
#' @keywords internal
column_dictionary <- function(file_path){
  # Get the index of each column name
  cNames <- get_header(large_file = file_path, colnames_only = TRUE)
  colDict <- setNames(seq(1,length(cNames)), cNames)
  return(colDict)
}


