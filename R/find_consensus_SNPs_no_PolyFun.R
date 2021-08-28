find_consensus_SNPs_no_PolyFun <- function(finemap_dat,
                                           verbose=T){
  printer("Identifying UCS and Consensus SNPs without PolyFun",v=verbose)
  newDF <- find_consensus_SNPs(finemap_dat,
                               exclude_methods = "POLYFUN_SUSIE",
                               sort_by_support = F)
  finemap_dat$Consensus_SNP_noPF <- newDF$Consensus_SNP
  finemap_dat$Support_noPF <- newDF$Support
  return(finemap_dat)
}
