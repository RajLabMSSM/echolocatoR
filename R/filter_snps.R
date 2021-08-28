#' Filter SNPs
#'
#' Filter SNps by MAF, window size, min/max position, maxmimum number of SNPs, or gene coordinates.
#' You can also explicitly remove certain variants.
#'
#' @family SNP filters
#' @keywords internal
#' @examples
#' data("BST1");
#' subset_DT <- filter_snps(subset_DT=BST1)
filter_snps <- function(subset_DT,
                        bp_distance=500000,
                        remove_variants=F,
                        min_POS=NA,
                        max_POS=NA,
                        max_snps=NULL,
                        min_MAF=NULL,
                        trim_gene_limits=F,
                        verbose=T){
  printer("FILTER:: Filtering by SNP features.",v=verbose)
  if(remove_variants!=F){
    printer("+ FILTER:: Removing specified variants:",paste(remove_variants, collapse=','), v=verbose)
    try({subset_DT <- subset(subset_DT, !(SNP %in% remove_variants) )})
  }
  # Trim subset according to annotations of where the gene's limit are
  if(trim_gene_limits!=F){
    subset_DT <- gene_trimmer(subset_DT=subset_DT,
                              gene=trim_gene_limits,
                              min_POS=min_POS,
                              max_POS=min_POS)
  }
  if(!is.null(max_snps)){
    subset_DT <- limit_SNPs(max_snps = max_snps,
                            subset_DT = subset_DT)
  }
  if(!is.null(min_MAF) & (!any(is.na(min_MAF))) ){
    if(any(min_MAF>0, na.rm = T) & "MAF" %in% colnames(subset_DT)){
      printer("+ FILTER:: Removing SNPs with MAF <",min_MAF,v=verbose)
      subset_DT <- subset(subset_DT, MAF>=min_MAF)
    }
  }
  # Limit range
  if(!is.null(bp_distance)){
    subset_DT <- assign_lead_SNP(new_DT = subset_DT,
                                 verbose = verbose)
    lead.snp <- subset(subset_DT, leadSNP)
    subset_DT <- subset(subset_DT,
                        POS >= lead.snp$POS - bp_distance &
                          POS <= lead.snp$POS + bp_distance)
  }
  if(!is.na(min_POS)){subset_DT <- subset(subset_DT, POS>=min_POS)}
  if(!is.na(max_POS)){subset_DT <- subset(subset_DT, POS<=max_POS)}
  printer("+ FILTER:: Post-filtered data:",paste(dim(subset_DT), collapse=" x "),v=verbose)
  return(subset_DT)
}



