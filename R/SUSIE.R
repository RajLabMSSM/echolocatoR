
# ***************** #
####    SUSIE    ####
# ***************** #

# * Notes on L parameter
# + L is the expected number of causal variants
# + Increasing L increases computational time
# + L=1: Gives a good amount of variation in PIP.
# + L=2: Warns "IBSS algorithm did not converge in 100 iterations!", but gives good variation in PIP.
# + L=3: Warns "IBSS algorithm did not converge in 100 iterations!". All PIPs 1s and 0s.
# + These results seem to be at least partially dependent on whether the ethnic composition of the LD matrix.
# * Notes on variance:
#   + If 'estimate_residual_variance' = TRUE _without_ providing 'var_y' _and_ L>1, susieR will throw error:
#   __"Estimating residual variance failed: the estimated value is negative"__
# + Running susieR with 'var_y = var(b)' provides _exactly_ the same results.
# * Statistical Terms:
#   + posterior inclusion probability (PIP)
# + coefficient estimate (Beta)
# + Effect allele frequency (EAF)
# + The I^2 statistic describes the percentage of variation across studies that seems not to be due to chance.



###------#### MAIN FUNCTION ###------####

#' Fine-map with SUSIE
#'
#' Sum of Single Effects (SuSiE): Iterative Bayesian Step-wise Selection
#'
#' \strong{Notes on convergence:}
#' \code{susieR} will often give the warning: \code{IBSS algorithm did not converge in 100 iterations!}.
#' This means the results might not necessarily be reliable.
#' There's several things you can try to avoid this:
#' \itemize{
#' \item{Increase \code{max_causal} (e.g. 5 => 10).}
#' \item{Increase \code{max_iter} (e.g. 100 => 1000), though this will take longer.}
#' \item{Decrease the locus window size, which will also speed up the algorithm but potentially miss causal variants far from the lead SNP.}
#' }
#'
#' @param max_causal The maximum number of non-zero effects (and thus causal variants).
#' @param rescale_priors If prior probabiltities are supplied,
#' rescale them from 0-1 (i.e. \code{rescaled_priors = priors / sum(priors)}).
#' @param plot_track_fit Record each iteration and make a GIF of the
#' fine-mapping algorithm learning the causal variants.
#' \strong{WARNING!:} Making this plot can take a long time if there's many iterations.
#' @inheritParams susieR::susie_suff_stat
#' @inheritParams finemap_pipeline
#' @source
#' \href{https://stephenslab.github.io/susieR/}{GitHub}
#' \href{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12388}{Publication}
#' @examples
#' data("BST1"); data("LD_matrix");
#' # LD_matrix <- readRDS("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
#' finemap_DT <- SUSIE(subset_DT=BST1, LD_matrix=LD_matrix)
SUSIE <- function(subset_DT,
                  LD_matrix,
                  dataset_type="GWAS",

                  # susieR default max_causal=L=10
                  max_causal=5,
                  # susieR default sample_size=n=<missing>
                  sample_size=NULL,
                  # susieR default prior_weights=NULL
                  prior_weights=NULL,
                  # susieR default PP_threshold=coverage=.95
                  PP_threshold=.95,
                  # PolyFun uses default `scaled_prior_variance=0.0001` (susieR default=0.2)
                  scaled_prior_variance=0.0001,# (previously 0.001)
                  # susieR default estimate_residual_variance=T
                  estimate_residual_variance=T,
                  # susieR default estimate_prior_variance=T
                  estimate_prior_variance=T,
                  # susieR default residual_variance=NULL
                  residual_variance=NULL,
                  # susieR default max_iter=100
                  max_iter=100,

                  rescale_priors=T,
                  plot_track_fit=F,
                  verbose=T){
  # if sample_size is NULL then SUSIE fails
  if(!"N" %in% names(subset_DT)){
      subset_DT <- get_sample_size(subset_DT)
  }
  sample_size <- max(subset_DT$N)

  # susie_vars <- get_var_y(subset_DT, dataset_type)

  printer("+ SUSIE:: max_causal =",max_causal, v=verbose)
  if(!is.null(prior_weights)){
    printer("+ SUSIE:: Utilizing prior_weights for",length(prior_weights),"SNPs.",v=verbose)
    if(rescale_priors){
      printer("+ SUSIE:: Rescaling priors from 0-1",v=verbose)
      prior_weights <- prior_weights / sum(prior_weights, na.rm = T)
    }
  }
  sub.out <- subset_common_snps(LD_matrix=LD_matrix,
                                finemap_dat=subset_DT,
                                fillNA = 0,
                                verbose = F)
  LD_matrix <- as.matrix(sub.out$LD)
  subset_DT <- sub.out$DT

  library(susieR)
  # SUSIE's authors "merge[d] susie_ss and susie_bhat to susie_suff_stat" in 11/2019.
  susie_version <- utils::packageVersion("susieR")
  if(length(find("susie_bhat"))==0){
    printer("+ SUSIE:: Using susie_suff_stat from susieR",paste0("v",susie_version),v=verbose)
    susie_func <- get("susie_suff_stat", asNamespace("susieR"))
  } else {
    printer("+ SUSIE:: Using susie_bhat from susieR",paste0("v",susie_version),v=verbose)
    susie_func <- get("susie_bhat", asNamespace("susieR"))
  }


  fitted_bhat <- susie_func(bhat = subset_DT$Effect,
                            shat = subset_DT$StdErr,
                            R = LD_matrix,
                            n = sample_size, # Number of samples/individuals in the dataset
                            L = max_causal, # maximum number of non-zero effects
                            ## NOTE: setting L == 1 has a strong tendency to simply return the SNP with the largest effect size.
                            scaled_prior_variance = scaled_prior_variance, # 0.1: Equates to "proportion of variance explained"
                            estimate_prior_variance = estimate_prior_variance, # default = FALSE
                            residual_variance = residual_variance,
                            # Raising max_iter can help susie converge
                            max_iter = max_iter,
                            # standardize = TRUE,
                            estimate_residual_variance = estimate_residual_variance, # TRUE

                            # IMPORTANT!! susieR uses the missing() function,
                            ## which means supplying var_y=NULL will give you errors!!!
                            ## When var_y is missing, it will be calculated automatically.
                            # var_y = var_y, # Variance of the phenotype (e.g. gene expression, or disease status)

                            # A p vector of prior probability that each element is non-zero
                            prior_weights = prior_weights,
                            coverage = PP_threshold,
                            track_fit = plot_track_fit,

                            verbose = F)

  if(plot_track_fit){
    try({susieR::susie_plot_iteration(fitted_bhat, n_causal, 'test_track_fit')})
  }
  printer("+ SUSIE:: Extracting Credible Sets...",v=verbose)
  ## Get PIP
  subset_DT$PP <- susieR::susie_get_pip(fitted_bhat)
  ## Get CS assignments
  CS_indices <- susieR::susie_get_cs(fitted_bhat)$cs
  susie_snps <- names(fitted_bhat$X_column_scale_factors)
  CS <- lapply(CS_indices, function(x){susie_snps[x]})
  CS_dict <- list()
  for(i in 1:length(CS)){
    for(s in CS[[i]]){
      CS_dict <- append(CS_dict, setNames(i,s))
    }
  }
  subset_DT$CS <- lapply(subset_DT$SNP, function(x){
    if(x %in% names(CS_dict) & subset(subset_DT, SNP==x)$PP>=PP_threshold){
      CS_dict[[x]]
      } else{0}}) %>% unlist()
  return(subset_DT)
}







###------#### SUSIER SUPPORT FUNCTIONS ###------####





# get_var_y <- function(subset_DT, dataset_type){
#   if(dataset_type=="GWAS" & "N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
#     printer("++ Computing phenotype variance...")
#     phenotype_variance <- var(c(rep(0, max(subset_DT$N_cases)),
#                                 rep(1, max(subset_DT$N_controls)))
#     )
#   } else if(dataset_type=="eQTL" & "Expression" %in% colnames(subset_DT)){
#     phenotype_variance <- var(subset_DT$Expression)
#   } else {
#     printer("++ Phenotype variance could not be calculated from this data.")
#     printer("    Estimating prior variance instead...")
#     phenotype_variance <- 1
#   }
#   return(list(phenotype_variance=phenotype_variance))
# }




