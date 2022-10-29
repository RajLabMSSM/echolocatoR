#' Argument list handler
#'
#' Infer whether the argument should be applied to all loci or matched by index.
#' @param arg Name of the argument to evaluate.
#' @param i Iterator index.
#' @param loci A vector of loci that are being iterated over.
#' @param use_names Whether to identify locus-specific value
#' by name instead of index.
#' @param error Throw an error (\code{TRUE})
#'  instead of a warning (\code{FALSE}).
#' @returns Selected argument value.
#' @family general
#' @keywords internal
arg_list_handler <- function(arg,
                             i,
                             loci,
                             use_names=FALSE,
                             error=TRUE){
  env <- parent.frame(n = 1)
  arg_vals <- get(arg, pos = env)
  #### Identify locus value by name ####
  if(isTRUE(use_names)){
    locus <- loci[[i]]
    if(!locus %in% names(arg_vals)){
      stp <- paste0("Value for locus ",shQuote(locus),
                    " not present in ",arg,".")
      if(isTRUE(error)){
        stop(stp)
      } else{
        wrn <- paste(stp,"Setting to NULL.")
        warning(wrn)
        return(NULL)
      }
    } else {
      return(arg_vals[[locus]])
    }
  }
  #### Identify locus value by index ####
  output <- if(length(arg_vals)>1){
    if(length(loci)!=length(arg_vals)){
      stp <- paste(arg,"must have a length equal to 1 or the number of loci")
      if(isTRUE(error)){
        stop(stp)
      } else {
        wrn <- paste(stp,"Setting to NULL.")
        warning(wrn)
        return(NULL)
      }
    } else {
      arg_vals[[i]]
    }
  } else{
    arg_vals
  }
  return(output)
}
