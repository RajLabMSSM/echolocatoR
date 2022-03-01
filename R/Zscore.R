
#' Compute Z-score
#'
#' Computes Z-score when you don't have the full vector,
#' but you have the necessary info about the full vector stored in \code{z.info}.
#'
#'These functions are necessary for \code{PAINTOR}.
#' @keywords internal
Zscore <- function(x, z.info){
  # Need to use the mean and standard deviation of the FULL dataset (i.e. all beta fomr the full summary stats file)
  sample.stdv <- z.info$sample.stdv
  sample.mean <- z.info$sample.mean
  z <- (x - sample.mean) / sample.stdv
  return(z)
}


#' Compute Z-score
#'
#' Computes Z-score when you have the full vector of values (not just a subset).
#'
#' These functions are necessary for \code{PAINTOR}.
#' @keywords internal
zscore <- function(vec){
  z <- scale(vec, center = TRUE, scale = TRUE)[,1]
  return(z)
}



