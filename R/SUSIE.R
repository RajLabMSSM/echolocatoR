
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

SUSIE <- function(subset_DT,
                  LD_matrix,
                  dataset_type="GWAS",
                  n_causal=5,
                  sample_size=NA,
                  var_y="estimate",
                  prior_weights=NULL,
                  PP_threshold=.95,
                  scaled_prior_variance=0.001,
                  estimate_residual_variance=F){
  # Sum of Single Effects (SuSiE): Iterative Bayesian Step-wise Selection
  # https://stephenslab.github.io/susieR/
  vars <- get_var_y(subset_DT, dataset_type)
  sample_size <- get_sample_size(subset_DT, sample_size)
  printer("+ SUSIE:: n_causal =",n_causal)
  if(!is.null(prior_weights)){
    printer("Utilizing prior_weights for",length(prior_weights),"SNPs.")
  }
  # library(susieR)
  # subset_DT <- data.table::fread("~/Desktop/Omer_example/Multi-finemap_results.txt")
  # LD_matrix <- readRDS("~/Desktop/Omer_example/UKB_LD.RDS")
  # sub.out <- subset_common_snps(LD_matrix, subset_DT)
  # LD_matrix <- sub.out$LD
  # subset_DT <- sub.out$DT
  # library(susieR)
  # SUSIE's authors "merge[d] susie_ss and susie_bhat to susie_suff_stat" in 11/2019.
  library(susieR)
  susie_func <- ifelse(length(find("susie_bhat"))==0,
                       susieR::susie_suff_stat,
                       susieR::susie_bhat)
  fitted_bhat <- susie_func(bhat = subset_DT$Effect,
                            shat = subset_DT$StdErr,
                            R = LD_matrix,
                            n = sample_size, # Number of samples/individuals in the dataset
                            L = n_causal, # we assume there are at *most* 'L' causal variables
                            ## NOTE: setting L == 1 has a strong tendency to simply return the SNP with the largest effect size.
                            scaled_prior_variance = scaled_prior_variance, # 0.1: Equates to "proportion of variance explained"
                            estimate_prior_variance = T, # default = FALSE
                            residual_variance = NULL,
                            # standardize = TRUE,
                            estimate_residual_variance = estimate_residual_variance, # TRUE
                            var_y = vars$phenotype_variance, # Variance of the phenotype (e.g. gene expression, or disease status)

                            # A p vector of prior probability that each element is non-zero
                            prior_weights = prior_weights,
                            coverage = PP_threshold,

                            verbose = F)

  # try({susieR::susie_plot_iteration(fitted_bhat, n_causal, 'test_track_fit')})
  printer("")
  printer("++ Extracting Credible Sets...")
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
  subset_DT$CS <- lapply(subset_DT$SNP, function(x){ if(x %in% names(CS_dict) & subset(subset_DT, SNP==x)$PP>=PP_threshold){ CS_dict[[x]] } else{0}}) %>% unlist()
  return(subset_DT)
}







###------#### SUSIER SUPPORT FUNCTIONS ###------####





get_var_y <- function(subset_DT, dataset_type){
  if(dataset_type=="GWAS" & "N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
    printer("++ Computing phenotype variance...")
    phenotype_variance <- var(c(rep(0, max(subset_DT$N_cases)),
                                rep(1, max(subset_DT$N_controls)))
    )
  } else if(dataset_type=="eQTL" & "Expression" %in% colnames(subset_DT)){
    phenotype_variance <- var(subset_DT$Expression)
  } else {
    printer("++ Phenotype variance could not be calculated from this data.")
    printer("    Estimating prior variance instead...")
    phenotype_variance <- NA
  }
  return(list(phenotype_variance=phenotype_variance))
}




