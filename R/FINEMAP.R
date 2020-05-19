


# ***************** #
####   FINEMAP   ####   
# ***************** # 
construct_FINEMAP_data <- function(results_path, 
                                   subset_DT,
                                   LD_matrix){
  ####### data.z ####### 
  printer("++ Formatting data.z file for FINEMAP")
  data.z <- subset_DT %>% dplyr::select(rsid=SNP, 
                                        chromosome=CHR, 
                                        position=POS,
                                        allele1=A1, 
                                        allele2=A2, 
                                        maf=MAF, 
                                        beta=Effect, # *required
                                        se=StdErr # *required
  )
  data.z$flip <- 0 # [optional] - flip==1, don't flip==0
  
  # !!! IMPORTANT !!!
  # Trim whitespaces
  ## Extra whitespace causes problems when you try to make space-delimited files
  # https://stackoverflow.com/questions/20760547/removing-whitespace-from-a-whole-data-frame-in-r
  cols_to_be_rectified <- names(data.z)[vapply(data.z, is.character, logical(1))]
  data.z <- data.z %>% mutate_at(vars(cols_to_be_rectified),
                                 funs(trimws) )
  
  
  ####### data.ld #######
  printer("++ Formatting LD Matrix for FINEMAP")
  ## The order of the SNPs in the dataset.ld must correspond to the order of variants in dataset.z.
  # load(file.path(results_path,"plink","LD_matrix.RData")) 
  
  # Filter 
  data.z <- subset(data.z, rsid %in% rownames(LD_matrix))
  ## This filters AND sorts LD_matrix by the order of rsids in data.z
  LD_filt <- LD_matrix[data.z$rsid, data.z$rsid] 
  
  # Write files
  ## MUST be space-delimited
  printer("+++ Writing FINEMAP z and ld files...")
  if( dim(data.z)[1]==dim(LD_filt)[1] ){
    # data.z
    data.z_path <- file.path(results_path,"FINEMAP","data.z")
    data.table::fwrite(data.z, data.z_path, sep = " ", nThread = 4) 
    # Sys.chmod(data.z_path, "777", use_umask = FALSE)
    # data.ld
    data.ld_path <- file.path(results_path,"FINEMAP","data.ld")
    data.table::fwrite(data.table:::as.data.table.matrix(LD_filt),
                       data.ld_path, sep=" ", quote = F, col.names = F, 
                       nThread = 4)
    # Sys.chmod(data.ld_path, "777", use_umask = FALSE)
  } else {warning("+ FINEMAP:: Summary statistics file (data.z) and LD matrix (data.ld) must contain the same number of SNPs.")}
}

construct_FINEMAP_master <- function(results_path,   
                                     n_samples,
                                     dataset_number=1,
                                     file.k=NA){ # [optional input]){
  printer("++ Constructing FINEMAP master file.")
  # For full list of parameters: http://www.christianbenner.com 
  header <- "z;ld;snp;config;cred;log;n_samples"
  # pathList <-  paste(c(file.z, file.ld, file.snp, file.config, file.log, n_samples), collapse=";")
  files <- c("data.z",  # [required input]
             "data.ld", # [required input]
             "data.snp", # [output]
             "data.config", # [optional output]
             "data.cred", # [optional output]
             "data.log"# [optional output]
  )
  if(!is.na(file.k)){ pathList <- append(pathList, file.k) } 
  paths_list <- paste(c(file.path("FINEMAP",files),n_samples), collapse = ";")  
  # Write master file
  dir.create(file.path(results_path, "FINEMAP"), recursive = T, showWarnings = F)
  data.table::fwrite(list(header,paths_list), file.path(results_path,"FINEMAP","master"), quote=F, sep="\n")
}


process_FINEMAP_results <- function(results_path, 
                                    subset_DT, 
                                    credset_thresh=.95,
                                    pvalue_thresh=.05,
                                    finemap_version="1.3"){
  # Import credible sets   
  if(finemap_version=="1.4"){
    # Annoying formatting differences between versions....
    # data.cred <- data.table::fread(file.path(results_path,"FINEMAP/data.cred"), 
    #                                 skip=4, na.strings = c("<NA>","NA"))
    # cred.cols <- grep("cred*", colnames(data.cred), value = T)
    # prob.cols <- grep("prob*", colnames(data.cred), value = T)
    # CS <- lapply(i:nrow(data.cred), function(i){
    #     rsids <- subset(data.cred, select=cred.cols)[i,]
    #   PP_vals <- subset(data.cred, select=prob.cols)[i,]  
    #   cred_sets <- data.table::data.table(SNP=unname( t(rsids)[,1] ), 
    #              PP=unname(t(PP_vals)[,1]), 
    #              Credible_Set=i) 
    #   return(cred_sets)
    # }) %>% data.table::rbindlist()
    # subset(CS, !is.na(SNP))
    printer("FINEMAP:: !!UNDER CONSTRUCTION!!")
    top_config <- data.table::fread(file.path(results_path,"FINEMAP/data.config")) 
    top_config <- subset(top_config, pvalue<pvalue_thresh)[1,]
    
  } else {
     top_config <- data.table::fread(file.path(results_path,"FINEMAP/data.config")) 
     top_config <- subset(top_config, prob>=credset_thresh)[1,]
  }
 
  Credible_Set <- strsplit(top_config$config, ",")[[1]]
  # Import snp-level results
  snp_level <- data.table::fread(file.path(results_path,"FINEMAP/data.snp"), sep=" ")
  # Merge with original data 
  subset_DT$Credible_Set <- ifelse(subset_DT$SNP %in% Credible_Set, 1, 0)
  subset_DT <- data.table:::merge.data.table(data.table::data.table(subset_DT), 
                                             data.table::data.table(subset(snp_level, select=c("rsid","prob")) ),
                                             by.x = "SNP", by.y="rsid")
  subset_DT <- subset_DT %>% dplyr::rename(PP=prob) %>% arrange(desc(Credible_Set))
  return(subset_DT)
}



FINEMAP <- function(subset_DT,
                    results_path,
                    LD_matrix,
                    FINEMAP_path="./echolocatoR/tools/FINEMAP/finemap_v1.3_MacOSX",
                    n_samples=NA,
                    n_causal=5,# Max number of allowed causal SNPs
                    model="cond", # cond (stepwise conditional search) vs. sss (stochastic shotgun search)
                    remove_tmps=T, 
                    server=F,
                    credset_thresh=.95, 
                    finemap_version="1.3.1"){ 
  # http://www.christianbenner.com
  
  # The stepwise conditional search starts with a causal configuration containing the 
  ## SNP with the lowest P-value alone and then iteratively adds to the causal configuration 
  ## the SNP given the highest posterior model probability until no further SNP yields
  ## a higher posterior model probability.
  if(is.na(n_samples) & "N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
    printer("+ FINEMAP:: Inferring sample size.")
    n_samples <- max(subset_DT$N_cases) + max(subset_DT$N_controls)
  }
  # Setup files
  construct_FINEMAP_master(results_path = results_path, n_samples = n_samples)
  construct_FINEMAP_data(results_path = results_path, 
                         subset_DT = subset_DT, 
                         LD_matrix = LD_matrix) 
  # Command line
  ## Example: 
  ## cmd <- paste(FINEMAP_path," --sss --in-files",file.path(dirname(FINEMAP_path),"example","master"), "--dataset 1 --n-causal-snps 5") 
  if(startsWith(getwd(),"/sc/")){server <- T}
  if(server){  
      FINEMAP_path <- paste0("ml finemap/",finemap_version," && finemap")
  } else {
    file.copy(from=FINEMAP_path, to=file.path(results_path))
    finemap_version <- "1.3"
    FINEMAP_path <- "./finemap_v1.3_MacOSX"
  }
 
  cmd <- paste("cd",results_path,"&&",
               FINEMAP_path,
               paste0("--",model),
               "--in-files",file.path("FINEMAP/master"),
               "--log",
               # --n-causal-snps  Option to set the maximum number of allowed causal SNPs
               # Default is 5 
               "--n-causal-snps",n_causal) 
  printer(cmd)
  system(cmd) 
  if(!server){ file.remove(file.path(results_path,"finemap_v1.3_MacOSX")) }
  # Process results
  finemap_DT <- process_FINEMAP_results(results_path, 
                                        subset_DT, 
                                        credset_thresh = credset_thresh, 
                                        finemap_version = finemap_version)
  
  # Remove tmp files
  if(remove_tmps){
    printer("+FINEMAP: Removing tmp files...") 
    tmp_files <- file.path(results_path,"FINEMAP",
                           c("data.ld",
                             "data.log_cond",
                             "data.log_sss", 
                             "data.z",
                             "master")) 
    suppressWarnings(file.remove(tmp_files))
  }
  return(finemap_DT)
}





