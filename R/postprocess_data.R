postprocess_data <- function(FINEMAP_DAT,
                             loci,
                             return_all,
                             verbose){
  if(return_all){
    messager("Returning results as nested list.",v=verbose)
    names(FINEMAP_DAT) <- loci
    FINEMAP_DAT[["merged_dat"]] <- data.table::rbindlist(
      lapply(FINEMAP_DAT, function(x){x$finemap_dat}),
      fill = TRUE,
      use.names = TRUE,
      idcol = "Locus"
    )
  } else {
    messager("Returning results as a merged data.table only.",v=verbose)
    FINEMAP_DAT <- data.table::rbindlist(FINEMAP_DAT, fill = TRUE)
    try({
      FINEMAP_DAT <- echodata::find_consensus_snps(
        dat = FINEMAP_DAT,
        credset_thresh = credset_thresh,
        consensus_thresh = consensus_thresh,
        verbose = verbose)
    })
  }
  return(FINEMAP_DAT)
}
