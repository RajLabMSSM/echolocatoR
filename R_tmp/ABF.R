# ***************** #
####     ABF     ####   
# ***************** #


ABF <- function(subset_DT, 
                PP_threshold=.95){
  # printer("Fine-mapping with ABF... \n")
  #data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset.txt") 
  finemap_DT <- coloc::finemap.abf(dataset = list(beta = subset_DT$Effect,
                                                  varbeta = subset_DT$StdErr^2, # MUST be squared
                                                  N = length(subset_DT$Effect),
                                                  s = subset_DT$proportion_cases,
                                                  snp = subset_DT$SNP,
                                                  MAF = subset_DT$MAF, 
                                                  type="cc")) 
  finemap_DT <- subset(finemap_DT, snp!="null") %>% 
    dplyr::rename(SNP=snp, PP=SNP.PP) %>%
    arrange(desc(PP))
  # Arbitarily assign SNPs with the top N probability as the top candidate causal SNPs
  # finemap_DT$Credible_Set <- c(rep(1,n_causal), rep(0,dim(finemap_DT)[1]-n_causal))
  # Any SNPs with a PP greater than the set threshold get included in the credible set
  finemap_DT$Credible_Set <- ifelse(finemap_DT$PP >= PP_threshold, 1, 0)
  finemap_DT <- data.table:::merge.data.table(x=data.table::data.table(subset_DT), 
                                              y=data.table::data.table(subset(finemap_DT, select=c("SNP","PP","Credible_Set")) ), 
                                              on="SNP")
  finemap_DT <- finemap_DT %>% arrange(desc(PP))
  return(finemap_DT)
} 

