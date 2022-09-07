#' Identify SNPs to condition on
#'
#' When running conditional analyses (e.g. \emph{GCTA-COJO}),
#' this functions automatically identifies SNP to condition on.
#' @inheritParams finemap_locus
#' @inheritParams finemap_loci
#' @family SNP filters
#'
#' @keywords internal
#' @returns Vector of SNPs
snps_to_condition <- function(conditioned_snps,
                              topSNPs,
                              loci){
  if(conditioned_snps=="auto"){
    lead_SNPs_DT <- subset(topSNPs, Locus %in% loci)
    # Reorder
    lead_SNPs_DT[order(factor(lead_SNPs_DT$Locus,levels= loci)),]
    return(lead_SNPs_DT$SNP)
  } else {return(conditioned_snps)}
}
