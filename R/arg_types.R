#' Argument types
#'
#' Categorize arguments of a given function by type.
#' @param fun Function to evaluate.
#' @param keep_types Argument types to keep in output.
#' @param verbose Print messages.
#' @keywords internal
#' @importFrom methods is
#' @return Named list of argument types.
arg_types <- function(fun,
                      keep_types=c("value","function",
                                   "NULL","no_default"),
                      verbose=FALSE){
  args <- formals(fun)
  types <- lapply(stats::setNames(names(args),
                                  names(args)),
                  function(n){
                    messager("Evaluating: ",n,v=verbose)
                    if (is.null(args[[n]])){
                      messager(n, "has a default value: NULL",v=verbose)
                      return("NULL")
                    } else if(any(!nzchar(args[[n]]) &
                                  is.name(args[[n]]))) {
                      messager(n, "has no default value",v=verbose)
                      return("no_default")
                    } else if(methods::is(args[[n]],"call")){
                      messager(n, "has a default value:",
                               shQuote(deparse(args[[n]])),
                               v=verbose)
                      return("function")
                    } else{
                      if (any(!nzchar(eval(args[[n]])))){
                        messager(n, "has a default value:",
                                 shQuote(deparse(eval(args[[n]]))),
                                 v=verbose)
                      }else{
                        messager(n, "has a default value:",
                                 shQuote(deparse(eval(args[[n]]))),
                                 v=verbose)
                      }
                      return("value")
                    }
                  })
  #### Filter by type ####
  types <- types[types %in% keep_types]
  return(types)
}
