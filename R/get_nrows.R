#' @family general
#' @keywords internal
get_nrows <- function(large_file){
  printer("+ Calculating the number of rows in",basename(large_file),"...")
  if(endsWith(large_file,".gz")){
    out <- system(paste("zcat",large_file,"| wc -l"), intern=T)
  } else {
    out <- system(paste("wc -l",large_file), intern=T)
  }
  file_nrows <- as.numeric(strsplit(out," ")[[1]][1])
  printer("++ File contains",file_nrows,"rows.")
  return(file_nrows)
}
