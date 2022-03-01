#' Identify SNPs to condition on
#'
#' When running conditional analyses (e.g. \emph{GCTA-COJO}),
#' this functions automatically identifies SNP to condition on.
#'
#' @family SNP filters
snps_to_condition <- function(conditioned_snps, top_SNPs, loci){
  if(conditioned_snps=="auto"){
    lead_SNPs_DT <- subset(top_SNPs, Locus %in% loci)
    # Reorder
    lead_SNPs_DT[order(factor(lead_SNPs_DT$Locus,levels= loci)),]
    return(lead_SNPs_DT$SNP)
  } else {return(conditioned_snps)}
}
