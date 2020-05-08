
# ***************** #
####    SUSIE    ####
# ***************** #


###------#### MAIN FUNCTION ###------####

SUSIE <- function(subset_DT,
                  LD_matrix,
                  dataset_type="GWAS",
                  n_causal=5,
                  sample_size=NA,
                  var_y="estimate",
                  prior_weights=NULL,
                  PP_threshold=.95){
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
  library(susieR)
  # SUSIE's authors "merge[d] susie_ss and susie_bhat to susie_suff_stat" in 11/2019.
  susie_func <- ifelse(length(find("susie_bhat"))==0,
                       susieR::susie_suff_stat, susieR::susie_bhat)
  fitted_bhat <- susie_func(bhat = subset_DT$Effect,
                            shat = subset_DT$StdErr,
                            R = LD_matrix,
                            n = sample_size, # Number of samples/individuals in the dataset
                            L = n_causal, # we assume there are at *most* 'L' causal variables
                            ## NOTE: setting L == 1 has a strong tendency to simply return the SNP with the largest effect size.
                            scaled_prior_variance = 0.001, # 0.1: Equates to "proportion of variance explained"
                            estimate_prior_variance = TRUE, # default = FALSE
                            residual_variance = NULL,
                            # standardize = TRUE,
                            estimate_residual_variance = F, # TRUE
                            var_y = vars$phenotype_variance, # Variance of the phenotype (e.g. gene expression, or disease status)

                            # A p vector of prior probability that each element is non-zero
                            prior_weights = prior_weights,
                            coverage = PP_threshold,

                            verbose = FALSE)

  # try({susieR::susie_plot_iteration(fitted_bhat, n_causal, 'test_track_fit')})
  printer("")
  printer("++ Extracting Credible Sets...")
  finemap_DT <- subset_DT
  ## Get PIP
  finemap_DT$PP <- susieR::susie_get_pip(fitted_bhat)
  ## Get CS assignments
  CS_indices <- susieR::susie_get_cs(fitted_bhat)$cs
  susie_snps <- names(fitted_bhat$X_column_scale_factors)
  Credible_Sets <- lapply(CS_indices, function(x){susie_snps[x]})
  CS_dict <- list()
  for(i in 1:length(Credible_Sets)){
    for(s in Credible_Sets[[i]]){
      CS_dict <- append(CS_dict, setNames(i,s))
    }
  }
  finemap_DT$Credible_Set <- lapply(finemap_DT$SNP, function(x){ if(x %in% names(CS_dict) & subset(finemap_DT, SNP==x)$PP>=PP_threshold){ CS_dict[[x]] } else{0}}) %>% unlist()


  # printer("++ Merging susieR results with original data")
  # finemap_DT <- data.table:::merge.data.table(x = data.table::data.table(subset_DT, key="SNP"),
  #                                             y = data.table::data.table(res_DT, key="SNP"), all.x = T)
  # finemap_DT <- finemap_DT %>% arrange(desc(PP))
  # Assign credible set #, or 0 to denote that it's not part of any credible set
  ## NOTE: if a SNP is part of more than one list, the top-ranked group to which is belong is used
  return(finemap_DT)
}







###------#### SUPPORT FUNCTIONS ###------####

printer <- function(..., v=T){if(v){print(paste(...))}}


get_sample_size <- function(subset_DT, sample_size=NA){
  if(is.na(sample_size)){
    if("N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
      sample_size <- max(subset_DT$N_cases) + max(subset_DT$N_controls)
    } else {
      sample_size <- 10000
      printer("++ No sample size variable detected...Defaulting to:",sample_size)
    }
  }
  return(sample_size)
}

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



check_credible <- function(subset_DT, fitted_bhat){
  ## *** IMPORTANT! ***: In case susieR cannot identify any credible set,
  # take the snps with the top 5 PIPs and provide a warning message. Interpret these snps with caution.

  # Credible_Set <- subset_DT[ as.numeric(strsplit( as.character(summary(fitted_bhat)$cs$variable) ,",")), ]$SNP
  printer("++",length(Credible_Set),"SNPs included in Credible Set")
  return(Credible_Set)
}
error_handling <- function(code) {
  tryCatch(code,
           error = function(c) {
             printer("--- ERROR ---")
             printer("****** Could NOT identify credible set. Default to SNPs with the top 5 PIPs ******")
             CS <- finemap_DT %>% arrange(desc(PIP))
             Credible_Set <- as.character(CS[1:5,]$SNP)
             return(Credible_Set)
           },
           warning = function(c) "Warning",
           message = function(c) "Message"
  )
}





