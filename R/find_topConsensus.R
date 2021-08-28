#' Find the top Consensus SNP
#'
#' Identify the \code{top_N} Consensus SNP(s) per Locus,
#' defined as the Consensus SNPs with the highest mean PP across all fine-mapping tools used.
#' @keywords internal
find_topConsensus <- function(dat,
                              top_N=1,
                              grouping_vars=c("Locus")){
  # CGet top consensus SNPs
  top.consensus <- (dat %>%
                      dplyr::group_by(.dots=grouping_vars) %>%
                      subset(Consensus_SNP) %>%
                      dplyr::arrange(-mean.PP) %>%
                      # dplyr::arrange(-IMPACT_score) %>%
                      # dplyr::mutate(topConsensus = ifelse(Consensus_SNP & mean.PP==max(mean.PP,na.rm = T),T,F)) %>%
                      dplyr::slice(top_N))$SNP %>% unique()
  dat$topConsensus <- dat$SNP %in% top.consensus
  return(dat)
}


