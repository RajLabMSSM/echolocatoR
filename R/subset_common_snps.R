#' Subset LD matrix and dataframe to only their shared SNPs
#'
#' Find the SNPs that are shared between an LD matrix and another data.frame with a `SNP` column.
#' Then remove any non-shared SNPs from both objects.
#'
#' @family SNP filters
#' @return data.frame
#' @keywords internal
#' @examples
#' data("BST1"); data('LD_matrix');
#' finemap_dat=BST1
#' out <- subset_common_snps(LD_matrix=LD_matrix, finemap_dat=BST1)
subset_common_snps <- function(LD_matrix,
                               finemap_dat,
                               fillNA=0,
                               verbose=F){
  printer("+ Subsetting LD matrix and finemap_dat to common SNPs...", v=verbose)
  # Remove duplicate SNPs
  LD_matrix <- data.frame(as.matrix(LD_matrix))
  LD_matrix <- LD.fill_NA(LD_matrix = LD_matrix,
                          fillNA = fillNA,
                          verbose = verbose)
  ld.snps <- unique(c(row.names(LD_matrix), colnames(LD_matrix)))

  # Remove duplicate SNPs
  finemap_dat <- finemap_dat[!base::duplicated(finemap_dat$SNP),]
  fm.snps <- finemap_dat$SNP
  common.snps <- base::intersect(ld.snps, fm.snps)
  if(length(common.snps)==0) stop("No overlapping RSIDs between LD_matrix and subset_DT")
  printer("+ LD_matrix =",length(ld.snps),"SNPs.", v=verbose)
  printer("+ finemap_dat =",length(fm.snps),"SNPs.", v=verbose)
  printer("+",length(common.snps),"SNPs in common.", v=verbose)
  # Subset/order LD matrix
  new_LD <- LD_matrix[common.snps, common.snps]

  # Subset/order finemap_dat
  finemap_dat <- data.frame(finemap_dat)
  row.names(finemap_dat) <- finemap_dat$SNP
  new_DT <- unique(data.table::as.data.table(finemap_dat[common.snps, ]))
  # Reassign the lead SNP if it's missing
  new_DT <- assign_lead_SNP(new_DT, verbose = verbose)
  # Check dimensions are correct
  if(nrow(new_DT)!=nrow(new_LD)){
    warning("+ LD_matrix and finemap_dat do NOT have the same number of SNPs.",v=verbose)
    warning("+ LD_matrix SNPs = ",nrow(new_LD),"; finemap_dat = ",nrow(finemap_dat), v=verbose)
  }
  return(list(LD=as.matrix(new_LD),
              DT=new_DT))
}



