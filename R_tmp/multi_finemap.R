#///////////////#
# Multi-finemap #
#///////////////#

find_consensus_SNPs <- function(finemap_DT,
                                verbose=T,
                                credset_thresh=.95,
                                consensus_thresh=2,
                                sort_by_support=T,
                                exclude_methods=NULL){
  printer("+ Identifying Consensus SNPs...",v=verbose)
  exclude_methods <- append(exclude_methods,"mean")
  # Find SNPs that are in the credible set for all fine-mapping tools
  CS_cols <- colnames(finemap_DT)[endsWith(colnames(finemap_DT),".Credible_Set")]
  CS_cols <- CS_cols[!(CS_cols %in% paste0(exclude_methods,".Credible_Set"))]
  if(consensus_thresh=="all"){consensus_thresh<-length(CS_cols)}
  printer("++ support_thresh =",consensus_thresh)
  # Get the number of tools supporting each SNP
  ## Make sure each CS is set to 1
  support_sub <- subset(finemap_DT, select = CS_cols) %>% data.frame()
  support_sub[sapply(support_sub, function(e){e>1})] <- 1
  finemap_DT$Support <- rowSums(support_sub, na.rm = T)
  finemap_DT$Consensus_SNP <- finemap_DT$Support >= consensus_thresh
  # Sort
  if(sort_by_support){
    finemap_DT <- finemap_DT %>% arrange(desc(Consensus_SNP), desc(Support))
  }

  # Calculate mean PP
  printer("+ Calculating mean Posterior Probability (mean.PP)...")
  PP.cols <- grep(".PP",colnames(finemap_DT), value = T)
  PP.cols <- PP.cols[!(PP.cols %in% paste0(exclude_methods,".PP"))]
  PP.sub <- subset(finemap_DT, select=c("SNP",PP.cols)) %>% data.frame()# %>% unique()
  PP.sub[is.na(PP.sub)] <- 0
  if(NCOL(PP.sub[,-1]) > 1){
    finemap_DT$mean.PP <- rowMeans(PP.sub[,-1], na.rm = T)
  } else{
    finemap_DT$mean.PP <- PP.sub[,-1]
  }
  finemap_DT$mean.Credible_Set <- ifelse(finemap_DT$mean.PP>=credset_thresh,1,0)

  # PP.sub %>% arrange(desc(mean.PP)) %>% head()
   printer("++",length(CS_cols),"fine-mapping methods used.")
   printer("++",dim(subset(finemap_DT,Support>0))[1],"Credible Set SNPs identified.")
   printer("++",dim(subset(finemap_DT,Consensus_SNP==T))[1],"Consensus SNPs identified.")
  return(finemap_DT)
}
# finemap_DT <- find_consensus_SNPs(finemap_DT)
# results_path <- "Data/GWAS/Kunkle_2019/PTK2B"
# subset_DT <- fread( file.path(results_path,"PTK2B_Kunkle_2019_subset.txt"))
# dataset_type <- "GWAS"
# load(file.path(results_path, "plink/LD_matrix.RData"))



# Fine-map using multiple fine-mapping tools
multi_finemap <- function(results_path,
                          fullSS_path,
                          finemap_method_list,
                          subset_DT,
                          dataset_type,
                          LD_matrix=NULL,
                          n_causal=5,
                          sample_size,
                          snp_col="SNP",
                          freq_col="Freq",
                          effect_col="Effect",
                          stderr_col="StdErr",
                          pval_col="P",
                          N_cases_col="N_cases",
                          N_controls_col="N_controls",
                          A1_col="A1",
                          A2_col="A2",
                          PAINTOR_QTL_datasets=NULL,
                          PP_threshold=.95){
  printer("++ Fine-mapping using multiple tools:", paste(finemap_method_list, collapse=", "))
  # finemap_method_list <- finemap_method_list[finemap_method_list!="COJO"]
  # Check overlap
  sub.out <- subset_common_snps(LD_matrix, subset_DT)
  LD_matrix <- sub.out$LD
  subset_DT <- sub.out$DT

  select_cols <- colnames(subset_DT)[!grepl(colnames(subset_DT),
                                            pattern = paste(c(finemap_method_list,"Support","Consensus_SNP"),collapse = "|"))]
  merged_DT <- subset(subset_DT, select = select_cols)

  for(i in 1:length(finemap_method_list)){
    m <- finemap_method_list[i];
    message("Multi-finemap:: ",m)
    DT <- null_DT <- data.table::data.table(SNP = merged_DT$SNP, Credible_Set = NA, Probability = NA);
    # DT <- tryCatch({
      # EXPRESSION
     try({
      DT <- finemap_method_handler(fullSS_path = fullSS_path,
                                   results_path = results_path,
                                   finemap_method = m,
                                   subset_DT = data.table::as.data.table(subset_DT),
                                   dataset_type = dataset_type,
                                   LD_matrix = LD_matrix,
                                   n_causal = n_causal,
                                   sample_size = sample_size,
                                   snp_col = snp_col,
                                   freq_col = freq_col,
                                   effect_col = effect_col,
                                   stderr_col = stderr_col,
                                   pval_col = pval_col,
                                   N_cases_col = N_cases_col,
                                   N_controls_col = N_controls_col,
                                   A1_col = A1_col,
                                   A2_col = A2_col,
                                   PAINTOR_QTL_datasets = PAINTOR_QTL_datasets,
                                   PP_threshold = PP_threshold)
     })
      # },
      # # WARNING
      # warning = function(w){printer("WARNING")},
      # # ERROR
      # error = function(e){
      #   message("SUSIE::Error:: Could not identify Credible Set.");
      #   return(null_DT)
      #   }
      # ) ## End tryCatch
      try({printer("++ Credible Set SNPS identified =",nrow(subset(DT,Credible_Set>0)) )})
      # Add results to method-specific columns
      printer("+++ Merging",m,"results with multi-finemap data.");
      value_var <- if(m=="COJO"){"Conditioned_Effect"}else{"PP"};
      DT_select <- subset(DT, select = c("SNP","Credible_Set",value_var) );
      # Rename columns according to method name
      cols <- colnames(DT_select);
      colnames(DT_select) <- c("SNP", paste(m, cols[2:length(cols)], sep="." ));

      # Merge new columns into DT
      merged_DT <- data.table:::merge.data.table(data.table::as.data.table(merged_DT),
                                                 data.table::as.data.table(DT_select),
                                                 by="SNP", all = T);

  }
  return(merged_DT)
}


####### Fine-mapping Handler #######
create_method_dir <- function(results_path,
                              finemap_method,
                              compress=T){
  method_dir <- file.path(results_path, finemap_method)
  # Make finemapping results folder
  dir.create(method_dir, recursive = T, showWarnings = F)
  # Return results file name
  dataset <- basename(dirname(results_path))
  locus <- basename(results_path)
  file_dir <- file.path(method_dir,
                        paste0(locus,"_",dataset,"_",finemap_method,".tsv",
                               ifelse(compress,".gz","")))
  return(file_dir)
}

save_finemap_results <- function(finemap_DT, file_dir){
  data.table::fwrite(finemap_DT, file_dir, sep = "\t", na = NA, quote = F, nThread = 4)
}


finemap_method_handler <- function(results_path,
                                   fullSS_path,
                                   finemap_method="SUSIE",
                                   subset_DT,
                                   dataset_type="GWAS",
                                   force_new_finemap=T,
                                   LD_matrix=NULL,
                                   n_causal=5,
                                   conditioned_snps,
                                   sample_size=NA,
                                   snp_col="SNP",
                                   freq_col="Freq",
                                   effect_col="Effect",
                                   stderr_col="StdErr",
                                   pval_col="P",
                                   N_cases_col="N_cases",
                                   N_controls_col="N_controls",
                                   A1_col="A1",
                                   A2_col="A2",
                                   PAINTOR_QTL_datasets=NULL,
                                   PP_threshold=.95){
  sub.out <- subset_common_snps(LD_matrix=LD_matrix,
                                finemap_DT=subset_DT)
  LD_matrix <- sub.out$LD
  subset_DT <- sub.out$DT

  message(paste("echolocatoR::",finemap_method))

  # INITIATE FINE-MAPPING
  if(finemap_method=="SUSIE"){
    # SUSIE
    finemap_DT <- SUSIE(subset_DT = subset_DT,
                        dataset_type = dataset_type,
                        LD_matrix = LD_matrix,
                        n_causal = n_causal,
                        # sample_size=N_cases+N_controls
                        sample_size = sample_size,
                        PP_threshold = PP_threshold)


  } else if(finemap_method=="POLYFUN_SUSIE"){
    # PolyFun+SUSIE
    finemap_DT <- POLYFUN_SUSIE(results_path = results_path,
                                finemap_DT = subset_DT,
                                LD_matrix = LD_matrix,
                                dataset_type = dataset_type,
                                n_causal = n_causal,
                                sample_size = sample_size,
                                polyfun_approach = "precomputed",#"non-parametric",
                                PP_threshold = PP_threshold)

  }else if(finemap_method=="ABF"){
    # coloc - finemap.abf
    finemap_DT <- ABF(subset_DT = subset_DT,
                      PP_threshold = PP_threshold)


  } else if(finemap_method=="FINEMAP"){
    # FINEMAP
    finemap_DT <- FINEMAP(subset_DT = subset_DT,
                          results_path = results_path,
                          LD_matrix = LD_matrix,
                          n_samples = sample_size,
                          n_causal = n_causal,
                          credset_thresh = PP_threshold)


  } else if("COJO" %in% finemap_method){
    #COJO
    conditioned_snps <- subset(subset_DT, leadSNP==T)$SNP
    finemap_DT <- COJO(subset_DT = subset_DT,
                       results_path = results_path,
                       fullSS_path = fullSS_path,
                       conditioned_snps = conditioned_snps,
                       conditional_analysis = T,
                       stepwise_procedure = F,

                       snp_col = snp_col,
                       freq_col = freq_col,
                       effect_col = effect_col,
                       stderr_col = stderr_col,
                       pval_col = pval_col,
                       N_cases_col = N_cases_col,
                       N_controls_col = N_controls_col,
                       A1_col = A1_col,
                       A2_col = A2_col)

  } else if("PAINTOR" %in% finemap_method) {
    finemap_DT <- PAINTOR(finemap_DT=subset_DT,
                           GWAS_datasets=ifelse(dataset_type=="GWAS",
                                                basename(dirname(results_path)),NULL),
                           QTL_datasets=NULL,
                           gene=basename(results_path),
                           n_causal=n_causal,
                           use_annotations=F,
                           PP_threshold=PP_threshold,
                           GWAS_populations="EUR",
                           LD_matrix=LD_matrix,
                           force_new_LD=F,
                           PP_threshold=PP_threshold)
  } else {
    stop("[::ERROR::] Enter valid finemap_method: 'SUSIE', 'ABF', 'FINEMAP', 'COJO', and 'PAINTOR' are currently available.")
  }
  return(finemap_DT)
}


finemap_handler <- function(results_path,
                            fullSS_path,
                            finemap_methods=c("SUSIE","FINEMAP"),
                            subset_DT,
                            dataset_type="GWAS",
                            force_new_finemap=T,
                            LD_matrix=NULL,
                            n_causal=5,
                            sample_size,
                            conditioned_snps,
                            snp_col="SNP",
                            freq_col="Freq",
                            effect_col="Effect",
                            stderr_col="StdErr",
                            pval_col="P",
                            N_cases_col="N_cases",
                            N_controls_col="N_controls",
                            A1_col="A1",
                            A2_col="A2",
                            PAINTOR_QTL_datasets=NULL,
                            PP_threshold=.95){
  message("-------- Step 4: Statistically Fine-map --------")
  start_FM <- Sys.time()
  set.seed(1)
  # First, check if there's more than one fin-mapping method given. If so, switch to multi-finemap function
  # if(length(finemap_methods)>1){
    ## Next, see if fine-mapping has previously been done (with multi-finemap)
    file_dir <- create_method_dir(results_path = results_path,
                                  finemap_method = "Multi-finemap",
                                  compress = T)
    old_file_dir <- file.path(dirname(file_dir),"Multi-finemap_results.txt")
    if(!file.exists(file_dir) & file.exists(old_file_dir)){file_dir <- old_file_dir }
    
    ### If so, import the previous results
    if(file.exists(file_dir) & force_new_finemap==F){
      printer("++ Previously multi-fine-mapped results identified. Importing...")
      finemap_DT <- data.table::fread(file_dir, sep="\t") %>% data.table::data.table()
    } else {
      ### If not, or if forcing new fine-mapping is set to TRUE, fine-map using multiple tools
      finemap_DT <- multi_finemap(results_path = results_path,
                                 fullSS_path = fullSS_path,
                                 finemap_method_list = finemap_methods,
                                 subset_DT = subset_DT,
                                 dataset_type = dataset_type,
                                 LD_matrix = LD_matrix,
                                 n_causal = n_causal,
                                 sample_size = sample_size,

                                 snp_col = snp_col,
                                 freq_col = freq_col,
                                 effect_col = effect_col,
                                 stderr_col = stderr_col,
                                 pval_col = pval_col,
                                 N_cases_col = N_cases_col,
                                 N_controls_col = N_controls_col,
                                 A1_col = A1_col,
                                 A2_col = A2_col,
                                 PAINTOR_QTL_datasets = PAINTOR_QTL_datasets,
                                 PP_threshold = PP_threshold)
      save_finemap_results(finemap_DT, file_dir)
    }
  # } else {
  #   # If only one tool is given, fine-map using only that tool
  #   ## Next, see if fine-mapping has previously been done (with this tool).
  #   file_dir <- create_method_dir(results_path, finemap_methods)
  #   ### If so, import the previous results
  #   if(file.exists(file_dir) & force_new_finemap==F){
  #     printer("\n +++ Previously fine-mapped", finemap_methods, "results identified. Importing...")
  #     finemap_DT <- data.table::fread(file_dir, sep="\t") %>% data.table::data.table()
  #   } else{
  #     ### If not, or if forcing new fine-mapping is set to TRUE, fine-map with that given tool
  #     finemap_DT <- finemap_method_handler(results_path = results_path,
  #                                  finemap_methods = finemap_methods,
  #                                  subset_DT = subset_DT,
  #                                  dataset_type = dataset_type,
  #                                  LD_matrix = LD_matrix,
  #                                  n_causal = n_causal,
  #                                  sample_size = sample_size)
  #     save_finemap_results(finemap_DT, file_dir)
  #   }
  # }
  end_FM <- Sys.time()
  printer("++ Fine-mapping with '", paste0(finemap_methods, collapse=", "),"' completed in ",round(end_FM-start_FM,2)," seconds.", sep="")
  return(finemap_DT)
}


#
# run_finemapping_tools <- function(X, Y, Z=NA){
#   # https://stephenslab.github.io/susie-paper/manuscript_results/pedagogical_example.html
#   data("N2finemapping")
#   attach(N2finemapping)
#   dat = N2finemapping
#   names(dat)
#
#   library(abind)
#   mm_regression = function(X, Y, Z=NULL) {
#     if (!is.null(Z)) {
#       Z = as.matrix(Z)
#     }
#     reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(susieR:::univariate_regression(X, Y[,i], Z)))
#     reg = do.call(abind, c(reg, list(along=0)))
#     # return array: out[1,,] is betahat, out[2,,] is shat
#     return(aperm(reg, c(3,2,1)))
#   }
#   sumstats = mm_regression(as.matrix(X), as.matrix(Y))
#   saveRDS(list(data=dat, sumstats=sumstats), 'N2.with_sumstats.rds')
#
#   printer("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH",
#       "Rscript ~/GIT/susieR/inst/code/finemap.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.FINEMAP\" args=\"--n-causal-max\ 2\" &> /dev/null")
#
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   Rscript ~/GIT/susieR/inst/code/caviar.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.CAVIAR\" args=\"-c\ 2\ -g\ 0.01\" &> /dev/null")
#
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   python ~/GIT/susieR/inst/code/dap-g.py N2.with_sumstats.rds N2finemapping.DAP -ld_control 0.20 --all &> /dev/null")
#
# }




