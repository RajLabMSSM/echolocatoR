set_tryCatch <- function(use_tryCatch=TRUE){
  if(isTRUE(use_tryCatch)){
    return(function(e){message(e);NULL})
  } else {
    return(function(e){stop(e)})
  }
}
