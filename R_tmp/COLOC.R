#     %%%%%%%%%%%%%%%%%      #
####### Colocalization ####### 
#     %%%%%%%%%%%%%%%%%     #

# Jansen et al. 2017: https://eqtl.onderzoek.io/index.php?page=gene_cis_details&gene=BST1

#-------------------------------------------------------
#-------------------------------------------------------
# ECHOLOCATOR DOCUMENTATION
# The input for coloc is a merged data.frame with the summary stats from both the GWAS and QTL, with one unique SNP per row. Those columns should be named as follows:
#   
#   ## Data.frame
#   # Shared columns
#   SNP: unique identifier for SNP that uses the same naming system for both the GWAS and QTL (e.g. rsid (preferred), or chr:position).
# # GWAS columns
# P: Nominal P-values from the GWAS. Technically you should be able to supply corrected P-values as well (e.g. FDR) since coloc will turn them into 
# Effect: Effect size (beta, or odds ratio) of each SNP in the GWAS
# StdErr: The standard error of each SNP in the GWAS (it should be non-squared since I square it in the coloc function).
# MAF: Minor allele frequency, preferably from the study itself but you can use a reference panel (e.g. 1000 Genomes) if necessary.
# # QTL columns
# QTL.P: Nominal P-values from the QTL.
# QTL.Effect: Effect size (beta, or odds ratio) of each SNP in the QTL
# QTL.StdErr: The standard error of each SNP in the QTL (it should be non-squared since I square it in the coloc function).
# QTL.MAF: Minor allele frequency, preferably from the study itself but you can use a reference panel (e.g. 1000 Genomes) if necessary, or even the MAF from the GWAS.
# 
# ## Single-value parameters
# sample_size: In addition to this data.frame, you need to supply the sample_size for both the GWAS and QTL, respectively. Rather than one value per SNP, these will be just two values (number of subjects in the GWAS and number of subjects in the QTL). Technically, sample_size is an optional parameter for coloc but I find it actually can have a large impact on the results, so supplying an accurate number for this is still important.
# proportion_cases: The proportion of subjects that were cases (e.g. PD). This is only necessary for GWAS, since QTL doesn't really have a case-control design.
#-------------------------------------------------------
#-------------------------------------------------------


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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

COLOC.construct_dataset <- function(subset_DT, 
                                    sample_size=NA, 
                                    proportion_cases=NA, # Doesn't allow actual 0, so use smallest number R allows
                                    MAF=NA,
                                    type="cc"){ 
  # Sample size
  sample_size <- get_sample_size(subset_DT, sample_size)
  # Proportion cases
  if(is.na(proportion_cases)){
    if("proportion_cases" %in% colnames(subset_DT)){
      printer("++ Extracting Proportion of Cases...")
      proportion_cases <- getmode(subset_DT$proportion_cases)
    }
  }
  # MAF
  if(length(MAF)==1){
    printer("++ Extracting MAF...")
    if(is.na(MAF) & "MAF" %in% colnames(subset_DT)){
      MAF <- subset_DT$MAF
    }  
  }
  # List form
  dataset <- list(pvalues = subset_DT$P, 
                  beta = subset_DT$Effect,
                  varbeta = subset_DT$StdErr^2, # MUST be squared
                  snp = subset_DT$SNP,
                  
                  N = sample_size, # [optional]
                  s = proportion_cases, # use overall proportions
                  MAF = MAF, # [required]
                  type = type)
  return(dataset)
}

 

COLOC <- function(data1,
                  data2,
                  data1_type="cc",# case-control
                  data2_type="quant",
                  shared_MAF=F, 
                  data1_proportion_cases=NA,
                  data2_proportion_cases=5e-324, 
                  data1_MAF="MAF",
                  data2_MAF=NA,
                  PP_threshold=0.8,
                  save_results=T,
                  results_path=NULL,
                  show_plot=T,
                  title1=NULL,
                  subtitle1=NULL,
                  title2=NULL,
                  subtitle2=NULL){
  printer("******** Step 5: COLOCALIZE ********")  
  # The Approximate Bayes Factor colocalisation analysis described in the next section 
  ## essentially works by fine mapping each trait under a single causal variant assumption 
  ##and then integrating over those two posterior distributions to calculate probabilities that 
  ## those variants are shared.
  # https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html
  # data1 <- data.table::fread(dataset1_path, stringsAsFactors = F, sep="\t") %>% data.frame()
  # data2 <- data.table::fread(dataset2_path, stringsAsFactors = F, sep="\t") %>% data.frame()
  common_SNPs <- intersect(data1$SNP, data2$SNP)
  data1 <- subset(data1, SNP %in% common_SNPs) %>% group_by(SNP) %>% slice(1) %>% data.frame()
  data2 <- subset(data2, SNP %in% common_SNPs) %>% group_by(SNP) %>% slice(1) %>% data.frame()

  dataset1 <- COLOC.construct_dataset(subset_DT = data1,
                                      sample_size = data2_sample_size,
                                      proportion_cases = data1_proportion_cases,
                                      type = data1_type, 
                                      MAF = data1_MAF)
  dataset2 <- COLOC.construct_dataset(subset_DT = data2, 
                                      sample_size = data2_sample_size,
                                      proportion_cases = data2_proportion_cases,
                                      type = data2_type, 
                                      MAF = data2_MAF)
  if(shared_MAF==1){
    data2$MAF <- data1$MAF
  } else if(shared_MAF==2){
    data1$MAF <- data2$MAF
  }
  ## NOTES: MESA and Fairfax: No sample size (SNP-level), proportion of cases, or freq/MAF info available?   
  # printer("\n\n")
  coloc.res <- coloc::coloc.abf(dataset1 = data1,
                                dataset2 = data2)
                                    # MAF = dataset1$MAF) 
 hypothesis_key <- setNames(
   c("Neither trait has a genetic association in the region.",
     "Only trait 1 has a genetic association in the region.",
     "Only trait 2 has a genetic association in the region.",
     "Both traits are associated, but with different causal variants (one in each dataset).",
     "Both traits are associated and share a single causal variant.") ,
   c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")) 
 # Report hypothess results
  printer("\n Hypothesis Results @ PP_threshold =",PP_threshold,":")
  true_hyp <-""
  for(h in names(hypothesis_key)){
    if(coloc.res$summary[h]>=PP_threshold){
      hyp <- hypothesis_key[h]
      printer("    ",h,"== TRUE: **",hyp )
      true_hyp <- paste0(names(hyp),": ", hyp)
    } else{
      printer("    ",h,"== FALSE: ")
    } 
  } 
   
  coloc_DT <- coloc.res$results 
  # Find the causal SNP that coloc.abf identified in each dataset via finemap.abf
  # DT1
  causal_DT1 <- coloc_DT %>% arrange(desc(lABF.df1))
  causal_DT1 <- causal_DT1$snp[1]
  # DT2
  causal_DT2 <- coloc_DT %>% arrange(desc(lABF.df2))
  causal_DT2 <- causal_DT2$snp[1]
  
  # Process results  
  coloc_DT$Colocalized <- ifelse(coloc_DT$SNP.PP.H4 >= PP_threshold, T, F)
  colocalized_snps <- subset(coloc_DT, Colocalized==T)$snp# subset(coloc_DT, Colocalized==1)$SNP
  coloc_datasets <- coloc_plot_data(coloc.res, data1, data2)
  subtitle2 <- paste0("Colocalized SNPs: ", paste(colocalized_snps,sep=", "))
  if((coloc.res$summary["PP.H3.abf"] + coloc.res$summary["PP.H4.abf"] >= PP_threshold) & 
     (coloc.res$summary["PP.H4.abf"]/coloc.res$summary["PP.H3.abf"] >= 2)){
    # "We called the signals colocalized when (coloc H3+H4 ≥ 0.8 and H4∕H3 ≥ 2)" -Yi et al. (2019)
    report <- paste("Datasets colocalized")  
  } else {
    report <- paste("Datasets NOT colocalized") 
  }   
  printer("\n++",report,"at PP.H3 + PP.H4 >=",PP_threshold," and PP.H3 / PP.H4 >= 2.","\n") 
 
  
  # Plot 
  if(show_plot){
    # title1 <- paste(strsplit(dataset1_path,"/")[[1]][4],strsplit(dataset1_path,"/")[[1]][3])
    # title2 <- paste(strsplit(dataset2_path,"/")[[1]][4],strsplit(dataset2_path,"/")[[1]][3]) 
    
    COLOC.plot(coloc_DT1 = coloc_datasets$coloc_DT1,
               coloc_DT2 = coloc_datasets$coloc_DT2, 
               title1 = title1,
               subtitle1 = subtitle1,
               title2 = title2, 
               subtitle2 = subtitle2,
               # SNP_list = c("rs34637584","rs76904798","rs117073808"),
               alt_color_SNPs = colocalized_snps, 
               show_plot = T)
  } 
  # Save
  if(save_results){
    coloc_path <- file.path(dirname(dataset1_path),"COLOC")
    dir.create(coloc_path, recursive = T, showWarnings = F)
    data.table::fwrite(coloc_DT, file.path(coloc_path,"COLOC_results.txt"), sep="\t")
    data.table::fwrite(coloc_DT, file.path(results_path, "COLOC/COLOC_raw.txt"), sep="\t")
  }
  return(coloc_DT)
} 



COLOC.plot_data <- function(coloc.res, 
                            data1, 
                            data2){
  coloc_DT <- coloc.res$results
  # Merge Dataset 1
  coloc_DT1 <- coloc_DT %>% dplyr::select(SNP="snp", 
                             V = "V.df1", 
                             Z= "z.df1", 
                             lABF = "lABF.df1", 
                             H4.Probability = "SNP.PP.H4")
  coloc_DT1 <- data.table:::merge.data.table(data1, coloc_DT1, on = "SNP", all = F)
  # Merge Dataset 2
  coloc_DT2 <- coloc_DT %>% dplyr::select(SNP="snp", 
                                          V = "V.df2", 
                                          Z= "z.df2", 
                                          lABF = "lABF.df2", 
                                          H4.Probability = "SNP.PP.H4")
  coloc_DT2 <- data.table:::merge.data.table(data2, coloc_DT2, on = "SNP", all = F) 
  return(list(coloc_DT1=coloc_DT1, coloc_DT2=coloc_DT2))
}

COLOC.plot <- function(gene, 
                       coloc_DT1, 
                       coloc_DT2, 
                       title1, 
                       title2, 
                       subtitle1,
                       subtitle2,
                       SNP_list=c(), 
                       alt_color_SNPs=c(), 
                       show_plot=T){ 
  mp1 <- manhattan_plot(subset_DT = coloc_DT1, 
                 gene = gene, 
                 SNP_list = SNP_list,
                 alt_color_SNPs = alt_color_SNPs,
                 title = title1,
                 subtitle = subtitle1, 
                 show_plot = F)
  mp2 <- manhattan_plot(subset_DT = coloc_DT2, 
                       gene = gene, 
                       SNP_list = SNP_list,
                       alt_color_SNPs = alt_color_SNPs,
                       title = title2,
                       subtitle = subtitle2, 
                       show_plot = F)
  cp <- cowplot::plot_grid(mp1, mp2, ncol = 1)
  if(show_plot){print(cp)}else{return(cp)}
}



COLOC.report_summary <- function(coloc.res, PP_threshold=.8){ 
  # MAF = dataset1$MAF) 
  hypothesis_key <- setNames(
    c("Neither trait has a genetic association in the region.",
      "Only trait 1 has a genetic association in the region.",
      "Only trait 2 has a genetic association in the region.",
      "Both traits are associated, but with different causal variants (one in each dataset).",
      "Both traits are associated and share a single causal variant.") ,
    c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")) 
  # Report hypothess results
  printer("Hypothesis Results @ PP_threshold =",PP_threshold,":")
  true_hyp <-""
  for(h in names(hypothesis_key)){
    if(!is.na(coloc.res$summary[h])){
      if(coloc.res$summary[h]>=PP_threshold){
        hyp <- hypothesis_key[h]
        printer("    ",h,"== TRUE: **",hyp )
        true_hyp <- paste0(names(hyp),": ", hyp)
      }
    } else{
      printer("    ",h,"== FALSE: ")
    } 
  } 
  
  # Save raw results   
  coloc_DT <- coloc.res$results
  # Process results  
  coloc_DT$Colocalized <- ifelse(coloc_DT$SNP.PP.H4 >= PP_threshold, T, F)
  colocalized_snps <- subset(coloc_DT, Colocalized==T)$snp# subset(coloc_DT, Colocalized==1)$SNP
  subtitle2 <- paste0("Colocalized SNPs: ", paste(colocalized_snps,sep=", "))
  if(!is.na(coloc.res$summary)["PP.H4.abf"] ){
    if((coloc.res$summary["PP.H3.abf"] + coloc.res$summary["PP.H4.abf"] >= PP_threshold) & 
       (coloc.res$summary["PP.H4.abf"]/coloc.res$summary["PP.H3.abf"] >= 2)){
      # "We called the signals colocalized when (coloc H3+H4 ≥ 0.8 and H4∕H3 ≥ 2)" -Yi et al. (2019)
      report <- paste("Datasets colocalized")  
    } else {report <- paste("Datasets NOT colocalized") }
  } else { report <- paste("Datasets NOT colocalized")}   
  printer("+ COLOC::",report,"at: PP.H3 + PP.H4 >=",PP_threshold," and PP.H3 / PP.H4 >= 2.") 
  return(coloc_DT)
}


# COLOC.plot_results <- function(coloc_DT){
#   # Find the causal SNP that coloc.abf identified in each dataset via finemap.abf
#   # DT1
#   causal_DT1 <- coloc_DT %>% arrange(desc(lABF.df1))
#   causal_DT1 <- causal_DT1$snp[1]
#   # DT2
#   causal_DT2 <- coloc_DT %>% arrange(desc(lABF.df2))
#   causal_DT2 <- causal_DT2$snp[1]
#   coloc_datasets <- coloc_plot_data(coloc.res, data1, data2)
#   # Plot 
#   title1 <- paste(gene,":",strsplit(dataset1_path,"/")[[1]][4],strsplit(dataset1_path,"/")[[1]][3])
#   title2 <- paste(gene,":",strsplit(dataset2_path,"/")[[1]][4],strsplit(dataset2_path,"/")[[1]][3]) 
#   
#   COLOC.plot(coloc_DT1 = coloc_datasets$coloc_DT1,
#              coloc_DT2 = coloc_datasets$coloc_DT2, 
#              title1 = title1,
#              subtitle1 = report,
#              title2 = title2, 
#              subtitle2 = subtitle2,
#              SNP_list = c("rs34637584","rs76904798","rs117073808"),
#              alt_color_SNPs = colocalized_snps, 
#              show_plot = T) 
# }


COLOC.PP4_plot <- function(COLOC_DT=NULL, coloc.results_path=NULL, PP_threshold=.8){
  if(is.null(COLOC_DT)){
    if(is.null(coloc.results_path)){coloc.results_path <- "./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/COLOC_results_noFlip-gwasEffect.txt"} 
    COLOC_DT <- data.table::fread(coloc.results_path)
  } 
  COLOC_DT$coloc <- (COLOC_DT$PP.H3.abf + COLOC_DT$PP.H4.abf >= PP_threshold) &  (COLOC_DT$PP.H4.abf/COLOC_DT$PP.H3.abf >= 2) 
  COLOC_DT <- COLOC_DT %>% dplyr::rename(QTL.Dataset=Dataset2)
  
  cp <- ggplot(subset(COLOC_DT, coloc==T), aes(x=QTL.Dataset, y=PP.H4.abf, fill=QTL.Dataset)) + 
    facet_grid(~Locus) + 
    geom_col(show.legend = F) + 
    coord_flip() + 
    theme_bw() + 
    scale_y_continuous(limits = c(0,1), breaks = c(0, 1)) +
    theme(rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))
  print(cp)
  ggsave(file.path(dirname(coloc.results_path),"COLOC_multiQTL.png"), 
         plot = cp, dpi=400, width = 14, height=5, bg="transparent")
}



COLOC.iterate_QTL <- function(GTEx_version="GTEx_V7",
                              dataset.gwas.name = "./Data/GWAS/Nalls23andMe_2019"){
  FM_all <- merge_finemapping_results(minimum_support = 0, 
                                      include_leadSNPs = T, 
                                      dataset = dataset.gwas.name, 
                                      xlsx_path = F)
  QTL_files <- list.files(path = "./Data/QTL", pattern = "*.finemap.txt.gz", recursive = T, full.names = T)
  QTL_files <- grep(paste(c("GTEx","Cardiogenics","MESA","psychENCODE"),collapse="|"), QTL_files, value = T)
  # QTL_files <- grep(paste(c("Fairfax"),collapse="|"), QTL_files, value = T)
 
  # qtl_file = QTL_files[28] # Need allele/MAF info: 1:4 (Brain_xQTL_Serve), 5:6 (Cardiogenics), 7:10 (Fairfax), 28:31 (psychENCODE) ##### (GTEx and MESA are good (tho MESA needs sample size))
  COLOC_DT <- lapply(QTL_files, function(qtl_file){
    dataset.qtl.name <- gsub(".finemap.txt.gz","",basename(qtl_file))
    printer("") 
    message("mergeQTL:: ", dataset.qtl.name)
    FM_merge <- mergeQTL.merge_handler(FM_all = FM_all, qtl_file = qtl_file)
    subset(FM_merge, SNP %in% c("rs7294619","rs76904798","rs11175620")) %>% createDT()
    ## WARNING: DON'T try to allele flip. Coloc does this for you.
    # FM_merge <- mergeQTL.flip_alleles(FM_merge) # NOTE!: Allele flipping makes a BIG difference.
   
    ## WARNING: Probably not valid to apply PP in this way since the model isn't intended to handle this. 
    ## Instead, use something like PAINTOR that takes SNP-level PPs.
    ## Create measures adjusted by Posterior Probabilities from fine-mapping
    # FM_merge$mean.PP <- rowMeans(subset(FM_merge, select=grep(".PP",colnames(FM_merge))))
    # FM_merge$Adjusted.Effect <- FM_merge$Effect * FM_merge$mean.PP
    
    coloc_dt <- lapply(unique(FM_merge$Gene), function(gene){
          printer("+ COLOC::",gene)  
          FM_gene <- subset(FM_merge, Gene==gene) 
          FM_gene <- FM_gene[!is.na(FM_gene$QTL.Effect),]
          default_results <- data.table::data.table(Locus=gene,
                                                    Dataset1=dataset.gwas.name, 
                                                    Dataset2=dataset.qtl.name, 
                                                    nSNPs_locus=nrow(FM_merge),
                                                    nSNPs_overlap=nrow(FM_gene),
                                                    PP.H0.abf=NA,
                                                    PP.H1.abf=NA,
                                                    PP.H2.abf=NA,
                                                    PP.H3.abf=NA,
                                                    PP.H4.abf=NA)
          if(nrow(FM_gene)>0){
            # GWAS
            dataset.gwas <- list(pvalues = FM_gene$P, 
                                 beta = FM_gene$Effect,
                                 varbeta = FM_gene$StdErr^2, # MUST be squared
                                 snp = FM_gene$SNP,
                                 
                                 N =  max(FM_gene$N_cases, na.rm = T) + max(FM_gene$N_controls, na.rm = T), # [optional]
                                 s = getmode(FM_gene$proportion_cases), # use overall proportions
                                 MAF = FM_gene$MAF, # [required]
                                 type = "cc")
            # QTL
            if(all(is.na(FM_gene$QTL.MAF))){FM_gene$QTL.MAF <- FM_gene$MAF; printer("+ COLOC:: No QTL.MAF given. Using GWAS.MAF instead.")} 
            FM_gene$QTL.MAF <- FM_gene$MAF ##############111111111111111
            dataset.qtl <- list(pvalues = FM_gene$QTL.P, 
                                beta = FM_gene$QTL.Effect,
                                varbeta = FM_gene$QTL.StdErr^2, # MUST be squared
                                snp = FM_gene$SNP,
                                
                                N = max(FM_gene$QTL.SampleSize), # [optional]
                                # s = NA, # Not used for quant studies
                                #sdY = ,# [optional] for a quantitative trait, the population standard deviation of the trait. if not given, it can be estimated from the vectors of varbeta and MAF
                                MAF = FM_gene$QTL.MAF, # [required]
                                type = "quant")
            
          # If regression coefficients and variances are available, it calculates Bayes factors for association at each SNP.           
          ## If only p values are available, it uses an approximation that depends on the SNP's MAF and ignores any uncertainty 
          ## in imputation. Regression coefficients should be used if available.
            coloc.res <- coloc::coloc.abf(dataset1 = dataset.gwas,
                                          dataset2 = dataset.qtl)
            COLOC.report_summary(coloc.res, PP_threshold = .8) 
            # dat <- data.table::data.table(coloc.res$results)
            dat <- data.table::data.table(t(coloc.res$summary))
            dat <- cbind(Locus=gene, 
                         Dataset1=dataset.gwas.name, 
                         Dataset2=dataset.qtl.name, 
                         nSNPs_locus=nrow(FM_merge),
                         dat) 
            dat <- dat %>% dplyr::rename(nSNPs_overlap=nsnps)
            return(dat)
          } else{
            printer("COLOC:: No overlapping SNPs.")
            return(default_results)
          }
    }) %>% data.table::rbindlist(fill=T) 
    return(coloc_dt)
  }) %>% data.table::rbindlist(fill=T)
 
  # Save
  if(save_results){
    coloc_path <- file.path(dataset.gwas.name,"_genome_wide","COLOC")
    dir.create(coloc_path, recursive = T, showWarnings = F)
    data.table::fwrite(COLOC_DT, file.path(coloc_path,"COLOC_results_noFlip-gwasEffect.txt"), sep="\t")
  } 
  COLOC.PP4_plot(COLOC_DT, PP_threshold = .95)
  return(COLOC_DT)
}

  



