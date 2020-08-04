# ***************** #
####     ABF     ####
# ***************** #

#' Fine-map with ABF
#'
#' Conduct statistical (non-functional) fine-mapping with approximate Bayes factor (ABF).
#'
#' @source
#' \itemize{
#' \item JB Maller et al., Bayesian refinement of association signals for 14 loci in 3 common diseases. Nature Genetics. 44, 1294–1301 (2012).
#' \item J Wakefield, A bayesian measure of the probability of false discovery in genetic epidemiology studies. American Journal of Human Genetics. 81, 208–227 (2007).
#' }
#' @keywords internal
#' @examples
#' data("BST1");
#' finemap_DT <- BST1
#' finemap_DT <- ABF(subset_DT=finemap_DT)
ABF <- function(subset_DT,
                PP_threshold=.95,
                case_control = TRUE){
  #data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset.txt")
  if(!"MAF" %in% colnames(subset_DT)) MAF <- NULL;
  if(!"N" %in% names(subset_DT) & is.null(sample_size)){
    ss_df <- get_sample_size(subset_DT)
    sample_size <- if(is.null(ss_df)) NULL else max(subset_DT$N, na.rm = T)
  }

  if( case_control == TRUE){
  finemap_dat <- coloc::finemap.abf(dataset = list(beta = subset_DT$Effect,
                                                  varbeta = subset_DT$StdErr^2, # MUST be squared
                                                  N = sample_size,
                                                  s = subset_DT$proportion_cases,
                                                  snp = subset_DT$SNP,
                                                  MAF = MAF,
                                                  type = "cc"))
  }else{
  finemap_dat <- coloc::finemap.abf(dataset = list(beta = subset_DT$Effect,
                                                  varbeta = subset_DT$StdErr^2, # MUST be squared
                                                  N = sample_size,
                                                  snp = subset_DT$SNP,
                                                  MAF = MAF,
                                                  type = "quant"))

  }
  finemap_dat <- subset(finemap_dat, snp!="null") %>%
    dplyr::rename(SNP=snp, PP=SNP.PP) %>%
    dplyr::arrange(desc(PP))
  # Arbitarily assign SNPs with the top N probability as the top candidate causal SNPs
  # finemap_dat$CS <- c(rep(1,n_causal), rep(0,dim(finemap_dat)[1]-n_causal))
  # Any SNPs with a PP greater than the set threshold get included in the credible set
  finemap_dat$CS <- ifelse(finemap_dat$PP >= PP_threshold, 1, 0)
  finemap_dat <- data.table:::merge.data.table(x=data.table::data.table(subset_DT),
                                              y=data.table::data.table(subset(finemap_dat, select=c("SNP","PP","CS")) ),
                                              on="SNP")
  finemap_dat <- finemap_dat %>% dplyr::arrange(desc(PP))
  return(finemap_dat)
}

