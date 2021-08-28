#' Infer sample size from summary stats
#'
#' @family general
#' @keywords internal
#' @examples
#' data("BST1")
#' BST1 <- finemap_DT
#' subset_DT <- get_sample_size(subset_DT = finemap_DT)
get_sample_size <- function(subset_DT,
                            sample_size=NULL,
                            effective_ss=T,
                            verbose=T){
  printer("+ Preparing sample_size (N) column",v=verbose)
  if(!"N" %in% colnames(subset_DT)){
    if(is.null(sample_size)){
      if("N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
        if(effective_ss){
          printer("++ Computing effective sample size.",v=verbose)
          subset_DT$N <- round(4.0 / (1.0/subset_DT$N_cases + 1.0/subset_DT$N_controls), digits = 0)
        }else {
          printer("++ Inferring sample size from max(N_cases) + max(N_controls):",sample_size,v=verbose)
          subset_DT$N <- subset_DT$N_cases + subset_DT$N_controls
        }
      } else {
        subset_DT$N <- NULL
        printer("++ `sample_size` not provided.",v=verbose)
      }
    } else {
      printer(paste0("++ Using `sample_size = ",sample_size,"` ",v=verbose) );
      subset_DT$N <- sample_size
    }
  } else {printer("+ `N` column already present in data.",v=verbose)}
  return(subset_DT)
}

