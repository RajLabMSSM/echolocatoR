#' Calculate effective sample size
#'
#' This is an often-cited approximation for the effective sample size in a case-control study.
#' i.e., this is the sample size required to identify an association with a
#'  quantitative trait with the same power as in your present study.
#'  (from email correpsondence with Omer Weissbrod)
#'
#' @param finemap_dat Preprocessed \emph{echolocatoR} locus subset file.
#' Requires the columns \strong{N_cases} and \strong{N_controls}.
#' @examples
#' data("BST1");
#' finemap_DT <- effective_sample_size(finemap_dat=BST1)
#' @keywords internal
effective_sample_size <- function(finemap_dat,
                                  sample_size=NULL,
                                  verbose=T){
  if(is.null(sample_size)){
    if(all(c("N_cases","N_controls") %in% colnames(finemap_dat))){
      finemap_dat$N <- round(4.0 / (1.0/finemap_dat$N_cases + 1.0/finemap_dat$N_controls), digits = 0 )
      printer("Calculating effective sample size (`N`) from `N_cases` and `N_controls`", v=verbose)
    }
  }else{
    finemap_dat$N <- sample_size
    printer(paste0("Using `sample_size = ",sample_size,"` ") )
  }
  return(finemap_dat)
}

