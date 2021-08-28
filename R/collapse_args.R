collapse_args <- function(args_list){
  # args_list <- list("--n-iterations"=5000,"--sss"=NULL)
  # OR
  # args_list <- "--n-iterations 5000 --sss"
  if(length(args_list)==0){
    return(NULL)
  }else {
    if(class(args_list)=="character"){
      return(args_list)
    } else {
      args_str <- lapply(names(args_list), function(x){
        paste(x, args_list[[x]])
      }) %>% paste(collapse=" ")
      return(args_str)
    }
  }
}
