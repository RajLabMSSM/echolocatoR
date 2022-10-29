#' Identify SNPs to condition on
#'
#' When running conditional analyses (e.g. \emph{GCTA-COJO}),
#' this functions automatically identifies SNP to condition on.
#' @inheritParams finemap_locus
#' @inheritParams finemap_loci
#' @family SNP filters
#'
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom methods is
#' @returns Vector of SNPs
snps_to_condition <- function(conditioned_snps,
                              topSNPs,
                              loci){
  Locus <- NULL;

  #### No conditioned SNPs ####
  if(is.null(conditioned_snps)) return(NULL)
  #### Automated conditioned SNPs ####
  if(all(tolower(conditioned_snps)=="auto")){
    lead_SNPs_DT <- subset(topSNPs, Locus %in% loci)
    # Reorder
    lead_SNPs_DT <- lead_SNPs_DT[
      order(factor(lead_SNPs_DT$Locus,levels= loci)),
    ]
    conditioned_snps <- lapply(
      X = stats::setNames(unique(lead_SNPs_DT$Locus),
                          unique(lead_SNPs_DT$Locus)),
      FUN = function(locus){
      d <- subset(lead_SNPs_DT, Locus==locus)
      if(nrow(d)>0) {
        return(unique(d$SNP))
      } else {
        return(NULL)
      }
    })
    return(conditioned_snps)

  #### User-specified conditioned SNPs ####
  } else {
    ## Get the conditioned SNPs for a specific locus
    conditioned_snps <- if(methods::is(conditioned_snps,"list")) {
      conditioned_snps2 <- lapply(stats::setNames(loci,
                                                  loci),
                                 function(locus){
        if(locus %in% names(conditioned_snps)){
          return(conditioned_snps[[locus]])
        } else {
          return(NULL)
        }
      })
    } else if(methods::is(conditioned_snps,"character")){
      ## Simply supply all conditioned SNPs to all loci
      conditioned_snps2 <- lapply(stats::setNames(loci,
                                                  loci),
                                  function(locus){
        return(unname(conditioned_snps))
      })
    } else {
      stp <- paste(
        "conditioned_snps must be one of:",
        paste("\n -",c("named list","character vector","NULL"),collapse = "")
      )
      stop(stp)
    }
    return(conditioned_snps2)
  }
}
