#' Bind stored files
#'
#' Rapidly read a list of files from storage and concatenate them by rows.
#'
#' @family general
#' @keywords internal
.rbind.file.list <- function(file.list,
                             verbose=T,
                             nCores=4){
  merged.dat <- parallel::mclapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x)
    return(dat)
  }, mc.cores = nCores) %>% data.table::rbindlist(fill=T)
  return(merged.dat)
}

