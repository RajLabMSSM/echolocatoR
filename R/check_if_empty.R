#' @family general
#' @keywords internal
check_if_empty <- function(file_path){
  rowCheck <- dim(data.table::fread(file_path, nrows = 2))[1]
  if(rowCheck==0){
    stop("No SNPs identified within the summary stats file that met your criterion. :o")
  } else {printer("+ Subset file looks good.")}
}

