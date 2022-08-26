#' @family directory functions
directory_info <- function(info_path=NULL,
                           dataset_name,
                           variable="fullSS.local"){
  Data_dirs <- data.table::fread(info_path)
  directory = subset(Data_dirs, Dataset==dataset_name, select=variable) %>%
    as.character()
  return(directory)
}
