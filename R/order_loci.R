#' Order loci
#'
#' Order loci by UCS size, or alphabetically.
#' @keywords internal
#' @examples
#' merged_DT <- echodata::get_Nalls2019_merged()
#' merged_DT2 <- echolocatoR:::order_loci(dat=merged_DT)
order_loci <- function(dat,
                       by_UCS_size=FALSE,
                       descending=TRUE,
                       verbose=FALSE){
  if(by_UCS_size){
    messager("+ Ordering loci by UCS size.",v=verbose)
    locus_order <- echoannot:::get_CS_counts(merged_DT = dat)
    dat$Locus <- factor(dat$Locus,
                        levels = locus_order$Locus,
                        ordered = TRUE)
  } else{
    messager("+ Ordering loci alphabetically.",v=verbose)
    if(descending){
      dat$Locus <- factor(dat$Locus,
                          levels = rev(sort(unique(dat$Locus))),
                          ordered = TRUE)
    } else {
      dat$Locus <- factor(dat$Locus,
                          levels = sort(unique(dat$Locus)),
                          ordered = TRUE)
    }
  }
  return(dat)
}
