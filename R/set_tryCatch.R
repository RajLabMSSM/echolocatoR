set_tryCatch <- function(use_tryCatch=TRUE,
                         locus=NULL){
  if(isTRUE(use_tryCatch)){
    return(function(e){
      loc_msg <- if(!is.null(locus)) paste0("[Locus: ",locus,"] ") else ""
      ## Extract the call stack for context
      calls <- sys.calls()
      ## Find the most informative call (skip tryCatch internals)
      fn_names <- vapply(calls, function(cl){
        fn <- tryCatch(as.character(cl[[1]]), error=function(x) "")
        if(length(fn) > 1) fn <- paste(fn, collapse="::")
        fn
      }, character(1))
      ## Filter to echoverse function calls
      echo_fns <- fn_names[grepl("^echo|^catalog|^download", fn_names)]
      origin <- if(length(echo_fns) > 0) {
        paste0(" in ", echo_fns[length(echo_fns)], "()")
      } else ""
      message("\n",loc_msg,"ERROR",origin,": ",conditionMessage(e),"\n")
      NULL
    })
  } else {
    return(function(e){stop(e)})
  }
}
