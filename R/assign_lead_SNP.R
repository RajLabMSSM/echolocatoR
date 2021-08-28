#' Assign a lead GWAS SNP to a locus
#'
#' If none of the SNPs in the data.frame have \code{leadSNP==T},
#' then sort by lowest p-value (and then highest Effect size) and assign the top SNP as the lead SNP.
#'
#' @param data.frame Fine-mapping results data.frame.
#' @return Fine-mapping results data.frame with new boolean \strong{leadSNP} column,
#'  indicating whether each SNPs is the lead GWAS SNP in that locus or not.
#' @keywords internal
assign_lead_SNP <- function(new_DT, verbose=T){
  if(sum(new_DT$leadSNP,na.rm = T)==0){
    printer("+ leadSNP missing. Assigning new one by min p-value.", v=verbose)
    top.snp <- head(arrange(new_DT, P, desc(Effect)))[1,]$SNP
    new_DT$leadSNP <- ifelse(new_DT$SNP==top.snp,T,F)
  }
  return(new_DT)
}


