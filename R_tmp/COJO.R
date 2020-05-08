


make_cojo_path <- function(results_path){
  cojo_path <- file.path(results_path, "COJO")
  dir.create(cojo_path, recursive = T, showWarnings = F)
  return(cojo_path)
}



COJO_conditional <- function(GCTA_path, 
                             results_path,  
                             conditioned_snps, 
                             min_MAF=0, 
                             diff_freq=0.1, 
                             excluded_path){
  cojo_path <- make_cojo_path(results_path)
  cojoFile_path <- file.path(cojo_path,"cojo-file.ma")
  # Conditioning on a single SNP
  # NOTE: Can't condition on lead LRRK2 SNP because it's freq is too different from the reference population.
  ## 2 SNP(s) have large difference of allele frequency among the GWAS summary data and the reference sample.
  # Create file of SNPs you want to condition on
  snp_list = paste(conditioned_snps, collapse="\n") 
  conditioned_path <- file.path(cojo_path,"cojo-cond.txt")
  data.table::fwrite(list(snp_list),conditioned_path, quote = F, sep=" ")
  # Command line
  cojo_cmd1 <- paste(GCTA_path,
                     " --bfile ",file.path(results_path,"plink/plink"),
                     " --cojo-file ",cojoFile_path,
                     " --maf ",min_MAF,
                     " --diff-freq ",diff_freq,
                     " --cojo-cond ",conditioned_path,  
                     " --exclude ",excluded_path,
                     " --out ",file.path(cojo_path,"cojo"), sep="")
  printer("\n + COJO conditional analysis -- Conditioning on:",paste(snp_list, collapse=", ") )
  system(cojo_cmd1)
}
# Conditional results
get_conditional_results <- function(cojo_path){  
  conditional_results <- data.table::fread(file.path(cojo_path,"cojo.cma.cojo")) 
  conditional_results <- subset(conditional_results, pC<=0.05) %>% 
    arrange(desc(bC))
    # rename(CHR=Chr, POS=bp, Effect=b, P=p, Probability=bC)
  # conditional_results$credible_set <- c(rep(T,5), rep(F,dim(conditional_results)[1]-5))  
  return(conditional_results)
}


COJO_stepwise <- function(subset_DT,
                          GCTA_path,
                          results_path, 
                          min_MAF, 
                          excluded_path){
  cojo_path <- make_cojo_path(results_path)
  cojo.ma <- file.path(cojo_path,"cojo-file.ma")
  # Stepwise selection procedure to identify independent SNPs
  cojo_cmd2 <- paste(GCTA_path,
                     " --bfile ",file.path(results_path,"plink/plink"),
                     " --cojo-file ",cojo.ma, 
                     " --maf ",min_MAF,
                     " --cojo-slct",
                     # --cojo-slct: Perform a stepwise model selection procedure 
                     ## to select independently associated SNPs. 
                     ## Results will be saved in a *.jma file 
                     ## with additional file *.jma.ldr showing 
                     ## the LD correlations between the SNPs.  
                     " --exclude ",excluded_path,
                     " --out ",file.path(cojo_path,"cojo"), sep="")
  printer("\n + COJO Stepwise Selection Procedure -- Identifying independent SNPs...")
  system(cojo_cmd2) 
}
# Stepwise results
get_stepwise_results <- function(cojo_path){
  ## "Perform a stepwise model selection procedure to select independently associated SNPs. 
  ## Results will be saved in a *.jma file with additional file *.jma.ldr showing the LD correlations between the SNPs."
  independent_snps <- data.table::fread(file.path(cojo_path,"cojo.jma.cojo")) 
  dim(independent_snps)
  independent_snps$LD_r2 <- independent_snps$LD_r^2
  # stepwise_LDmatrix <-  data.table::fread(file.path(cojo_path,"cojo.ldr.cojo"), key = "SNP")
  return(independent_snps)
}

process_COJO_results <- function(subset_DT,
                                 results_path, 
                                 freq_cutoff=0.1){
  printer("+ Processing COJO results...")
  cojo_path <- make_cojo_path(results_path) 
  # Stepwise results 
  ## FILTER BY FREQ
  stepwise_results <- get_stepwise_results(cojo_path )%>% 
    subset(freq_geno > freq_cutoff)
  independent_SNPs <- stepwise_results$SNP
  # Conditional results
  cond_snps <- data.table::fread(file.path(cojo_path,"cojo-cond.txt"), header=F)$V1
  conditional_results <- get_conditional_results(cojo_path) %>% 
    dplyr::select(SNP, Conditioned_Effect = bC)
  
  colnames(conditional_results)[2] <- paste0(colnames(conditional_results)[2])
  # Merge with original data
  cojo_DT <- data.table:::merge.data.table(data.table::data.table(subset_DT, key = "SNP"),
                                           data.table::data.table(conditional_results, key = "SNP"),
                                              by = "SNP", all = T)
  cojo_DT$Credible_Set <- ifelse(cojo_DT$SNP %in% independent_SNPs,1,0)
  return(cojo_DT)
}


COJO <- function(subset_DT,
                 fullSS_path, 
                 results_path,
                 conditioned_snps,
                 excluded_snps="",
                 min_MAF=0, 
                 GCTA_path="echolocatoR/tools/gcta_1.92.1beta5_mac/bin/gcta64",
                 bfiles="plink_tmp/plink",
                 stepwise_procedure=T,
                 conditional_analysis=T, 
                 snp_col="SNP",
                 freq_col="Freq",
                 effect_col="Effect",
                 stderr_col="StdErr",
                 pval_col="P",
                 N_cases_col="N_cases",
                 N_controls_col="N_controls",
                 A1_col="A1",
                 A2_col="A2",
                 full_genome=F
                 ){ 
  ###########################################################
  # DOCUMENTATION: http://cnsgenomics.com/software/gcta/#COJO
  ###########################################################
  # 'Columns are SNP, the effect allele, the other allele, frequency of the effect allele, 
  ## effect size, standard error, p-value and sample size. 
  ## The headers are not keywords and will be omitted by the program. 
  ## Important: "A1" needs to be the effect allele
  ## with "A2" being the other allele and "freq" should be the frequency of "A1".'
  
  # Note: 1) For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error. 
  ## 2) Please always input the summary statistics of all SNPs even if your analysis only focuses on a subset of SNPs 
  ## because the program needs the summary data of all SNPs to calculate the phenotypic variance. 
  ## You can use one of the --extract options (Data management) to limit the COJO analysis in a certain genomic region.
  
  # Make path 
  cojo_path <- make_cojo_path(results_path) 
   
  ## Import SS subset
  # awk -F '\t' 'NR==1{print "SNP","A1","A2","freq","b","se","p","N"} NR>1{if($1==12 && $2>=40114434 && $2<=40935639){print $3, $4, $5, $6, $7, $8, $9, $10}}' nallsEtAl2019_allSamples_allVariants.mod.txt >  LRRK2_COJO.ma
  
  # Create cojo .ma file 
  ## NOTE: cojo-file.ma Must be a SPACE-separated file 
  if(full_genome){
    cojo.ma <- data.table::fread(fullSS_path) %>%
      dplyr::rename(N_cases = N_cases_col, N_controls = N_controls_col) %>%
      dplyr::mutate(N = N_cases + N_controls) %>% 
      dplyr::select(SNP=snp_col, 
                    A1, A2, 
                    freq=freq_col, 
                    b=effect_col, 
                    se=stderr_col, 
                    p=pval_col, 
                    N) 
    # Create genome-wide dir
    genome_dir <- file.path(dirname(results_path),"_genome_wide","COJO")
    dir.create(genome_dir, showWarnings = F, recursive = T)
    cojo_path <- genome_dir
    
# <<<<<<< HEAD:echolocatoR/R/Finemapping/COJO.R
#     # Convert all full vcfs to merged plink
#     base_dir <- "./Data/Reference/1000_Genomes/Phase1"
#     vcf_list <- list.files(base_dir, pattern = "*.vcf")
#     vcf_list <- fread(paste0(dirname(base_dir),"/1KG-P3_vcfs.txt"), header = F)$V1
#     
#     ## Convert
#     chromo_list <- c()
#     for(vcf in vcf_list){
#       chromo <- strsplit(basename(vcf),"[.]")[[1]][[2]]
#       cmd <- paste(plink_file(), 
#             " --vcf",vcf, 
#             "--maf 0.0", 
#             "--recode",
#             "--out",paste0(base_dir, chromo) 
#             )
#       chromo_list <- c(chromo_list, chromo)
#       # system(cmd)
#     }
#     
#     ## Merge
#     allfiles <- data.table::data.table(PED =  paste0(chromo_list,".ped"), 
#                                        MAP =  paste0(chromo_list,".map")) 
#     data.table::fwrite(allfiles,file.path(dirname(base_dir),"allfiles.txt"), sep="\t", col.names = F)
#     paste(plink_file(),
#           "--merge-list","allfiles.txt",
#           "--make-bed",
#           "--out MERGED")
#     
# =======
#     # Convert full vcf to plink
#     paste(plink_file(), " --vcf",vcf,
#     "myvcf.vcf --maf 0.05 --recode --out myplink")
# >>>>>>> 1e2aecb9b38f6c049a9c6f1d9baed0f0d268e0b4:echolocatoR/R/conditional.R
  } else {
    # Use subset of summary stats (not for the stepwise conditional procedure)
    cojo.ma <- subset_DT %>% dplyr::rename(N_cases = N_cases_col, 
                                           N_controls = N_controls_col, 
                                           A1 = A1_col,
                                           A2 = A2_col) %>%
      dplyr::mutate(N = N_cases + N_controls) %>% 
      dplyr::select(SNP, 
                    A1,
                    A2, 
                    freq="Freq", 
                    b="Effect", 
                    se="StdErr", 
                    p="P", 
                    N) 
  } 
 
  cojoFile_path <- file.path(cojo_path,"cojo-file.ma")
  data.table::fwrite(cojo.ma, cojoFile_path, sep=" ") 
  
  # Create of SNPs to exclude from analysis
  excluded_path <- file.path(cojo_path,"excluded_snps.txt")
  data.table::fwrite(list(excluded_snps),excluded_path , quote = F)
  
  ## Run COJO  
  if(stepwise_procedure){ 
    COJO_stepwise(GCTA_path = GCTA_path, 
                  results_path = results_path,  
                  min_MAF = min_MAF, 
                  excluded_path = excluded_path)
  }
  if(conditional_analysis){
    COJO_conditional(GCTA_path = GCTA_path, 
                     results_path = results_path, 
                     conditioned_snps = conditioned_snps, 
                     min_MAF = min_MAF, 
                     excluded_path = excluded_path)
  }
  
  cojo_DT <- process_COJO_results(subset_DT = subset_DT,
                                  results_path = results_path,
                                  freq_cutoff = 0.1)
  return(cojo_DT)
} 


