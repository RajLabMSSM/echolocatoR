#' Order loci by UCS size, or alphabetically
#'
#' @examples
#' data("merged_DT");
#' ... by UCS size ...
#' merged_dat <- order_loci(dat=merged_DT, merged_dat=merged_DT, descending=F)
#' ... alphabetically ...
order_loci <- function(dat,
                       merged_dat,
                       by_UCS_size=F,
                       descending=T,
                       verbose=F){
  if(by_UCS_size){
    printer("+ Ordering loci by UCS size.",v=verbose)
    locus_order <- SUMMARISE.get_CS_counts(merged_dat)
    dat$Locus <- factor(dat$Locus,  levels = locus_order$Locus, ordered = T)
  } else{
    printer("+ Ordering loci alphabetically.",v=verbose)
    if(descending){
      dat$Locus <- factor(dat$Locus,  levels = rev(sort(unique(merged_dat$Locus))), ordered = T)
    } else {
      dat$Locus <- factor(dat$Locus,  levels = sort(unique(merged_dat$Locus)), ordered = T)
    }
  }
  return(dat)
}

