#' Replace items within DT object
#'
#' Annoyingly, there is no native function to do simple find-and-replace in the `DT` library.
#'
#' @family general
#' @return data.frame
#' @keywords internal
dt_replace <- function(DT, target, replacement){
  for(col in names(DT)) set(DT, i=which(DT[[col]]==target), j=col, value=replacement)
  return(DT)
}
