#' Check required arguments
#'
#' Check whether any required arguments to a function are missing.
#' An argument is considered required if there is no default.
#' @inheritParams check_deprecated
#' @keywords internal
#' @returns Null
check_required <- function(fun=echolocatoR::finemap_loci,
                           args=match.call(call = sys.call(sys.parent(2)),
                                           expand.dots = FALSE),
                           return_no_default=FALSE){

  no_default <- arg_types(fun = fun,
                          keep_types = "no_default")
  missing_no_default <- no_default[!no_default %in% names(args)[-1]]
  if(length(missing_no_default)>0){
    stp <- paste("Missing",length(missing_no_default),"required argument(s):",
                 paste0("\n - ",names(missing_no_default),collapse = "")
                 )
    stop(stp)
  }
  if(isTRUE(return_no_default)) return(no_default)
}
