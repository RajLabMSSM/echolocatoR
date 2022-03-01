#' @family general
#' @keywords internal
arg_list_handler <- function(arg, i){
  output <- if(length(arg)>1){arg[i]}else{arg}
  return(output)
}
