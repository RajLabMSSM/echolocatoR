#' @family general
#' @keywords internal
get_nrows <- function(large_file){
  messager("+ Calculating the number of rows in",basename(large_file),"...")
  if(endsWith(large_file,".gz")){
    out <- system(paste("zcat",large_file,"| wc -l"), intern=TRUE)
  } else {
    out <- system(paste("wc -l",large_file), intern=TRUE)
  }
  file_nrows <- as.numeric(strsplit(out," ")[[1]][1])
  messager("++ File contains",file_nrows,"rows.")
  return(file_nrows)
}
