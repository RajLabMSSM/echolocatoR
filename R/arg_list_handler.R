#' Argument list handler
#'
#' Infer whether the argument should be applied to all loci or matched by index.
#' @param arg Name of the argument to evaluate.
#' @param i Iterator index.
#' @param loci A vector of loci that are being iterated over.
#' @returns Selected argument value.
#' @family general
#' @keywords internal
arg_list_handler <- function(arg,
                             i,
                             loci){
  env <- parent.frame(n = 1)
  arg_vals <- get(arg, pos = env)
  output <- if(length(arg_vals)>1){
    if(length(loci)!=length(arg_vals)){
      stp <- paste(arg,"must have a length equal to 1 or the number of loci")
      stop(stp)
    } else {
      arg_vals[[i]]
    }
  } else{
    arg_vals
  }
  return(output)
}
