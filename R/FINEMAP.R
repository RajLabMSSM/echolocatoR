


# ***************** #
####   FINEMAP   ####
# ***************** #


#' Prepare input files for \code{FINEMAP}
#'
#' Creates and saves 1) the summary stats file, and 2) the LD matrix.
#' @family FINEMAP
#' @keywords internal
#' @source
#' \url{http://www.christianbenner.com}
#' @examples
#' data("locus_dir"); data("BST1"); data("BST1_LD_matrix");
#' finemap_DT <- BST1
#' dir.create(file.path(locus_dir,"FINEMAP"), showWarnings = FALSE, recursive = TRUE)
#' out <- subset_common_snps(LD_matrix=LD_matrix, finemap_DT=finemap_DT)
#' LD_matrix <- out$LD
#' finemap_DT <- out$DT
#' dat_paths <- FINEMAP.construct_data(locus_dir=locus_dir, subset_DT=finemap_DT, LD_matrix=LD_matrix)
FINEMAP.construct_data <- function(locus_dir,
                                   subset_DT,
                                   LD_matrix,
                                   nThread=4){
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
  data.z <- data.z %>% mutate_at(.vars = vars(cols_to_be_rectified),
                                 .funs = trimws )

  ####### data.ld #######
  printer("++ Formatting LD Matrix for FINEMAP")
  ## The order of the SNPs in the dataset.ld must correspond to the order of variants in dataset.z.
  # load(file.path(locus_dir,"plink","LD_matrix.RData"))

  # Filter
  data.z <- subset(data.z, rsid %in% rownames(LD_matrix))
  ## This filters AND sorts LD_matrix by the order of rsids in data.z
  LD_filt <- LD_matrix[data.z$rsid, data.z$rsid]

  # Write files
  ## MUST be space-delimited
  printer("+++ Writing FINEMAP z and ld files...")
  if( dim(data.z)[1]==dim(LD_filt)[1] ){
    # data.z
    data.z_path <- file.path(locus_dir,"FINEMAP","data.z")
    data.table::fwrite(data.z, data.z_path, sep = " ",
                       nThread = nThread)
    # Sys.chmod(data.z_path, "777", use_umask = FALSE)
    # data.ld
    data.ld_path <- file.path(locus_dir,"FINEMAP","data.ld")
    data.table::fwrite(data.table:::as.data.table.matrix(LD_filt),
                       data.ld_path, sep=" ", quote = F, col.names = F,
                       nThread = nThread)
    # Sys.chmod(data.ld_path, "777", use_umask = FALSE)
  } else {warning("+ FINEMAP:: Summary statistics file (data.z) and LD matrix (data.ld) must contain the same number of SNPs.")}
  return(c("Zscore_path"=data.z_path,
           "LD_path"=data.ld_path))
}




#' Construct the \code{FINAMAP} master file
#'
#' Creates and saves the master file
#' which tells \code{FINEMAP} where to find each input file.
#' @family FINEMAP
#' @keywords internal
#' @source
#' \url{http://www.christianbenner.com}
#' @examples
#' data("locus_dir");
#' master_path <- FINEMAP.construct_master(locus_dir=locus_dir, n_samples=25000)
FINEMAP.construct_master <- function(locus_dir,
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
  dir.create(file.path(locus_dir, "FINEMAP"), recursive = T, showWarnings = F)
  master_path <- file.path(locus_dir,"FINEMAP","master")
  data.table::fwrite(list(header,paths_list), master_path, quote=F, sep="\n")
  return(master_path)
}




#' Post-processing of \code{FINEMAP} results
#'
#' @family FINEMAP
#' @keywords internal
#' @source
#' \url{http://www.christianbenner.com}
#' @examples
#' data("locus_dir"); data("BST1");
#' finemap_DT <- BST1
#' subset_DT <-FINEMAP.process_results(locus_dir=locus_dir, subset_DT=finemap_DT)
FINEMAP.process_results <- function(locus_dir,
                                    subset_DT,
                                    credset_thresh=.95,
                                    pvalue_thresh=.05,
                                    finemap_version="1.3"){
  # Import credible sets
  if(finemap_version=="1.4"){
    # Annoying formatting differences between versions....
    # data.cred <- data.table::fread(file.path(locus_dir,"FINEMAP/data.cred"),
    #                                 skip=4, na.strings = c("<NA>","NA"))
    # cred.cols <- grep("cred*", colnames(data.cred), value = T)
    # prob.cols <- grep("prob*", colnames(data.cred), value = T)
    # CS <- lapply(i:nrow(data.cred), function(i){
    #     rsids <- subset(data.cred, select=cred.cols)[i,]
    #   PP_vals <- subset(data.cred, select=prob.cols)[i,]
    #   cred_sets <- data.table::data.table(SNP=unname( t(rsids)[,1] ),
    #              PP=unname(t(PP_vals)[,1]),
    #              CS=i)
    #   return(cred_sets)
    # }) %>% data.table::rbindlist()
    # subset(CS, !is.na(SNP))
    printer("FINEMAP:: !!UNDER CONSTRUCTION!!")
    top_config <- data.table::fread(file.path(locus_dir,"FINEMAP/data.config"),nThread=nThread)
    top_config <- subset(top_config, pvalue<pvalue_thresh)[1,]

  } else {
    ## Configuration file
    ### Gives all model results for all the configurations tested
    ### (regardless of whether they're over the 95% probability threshold)
    # top_config <- data.table::fread(file.path(locus_dir,"FINEMAP/data.config"), nThread=nThread)
    # top_config <- subset(top_config, prob>=credset_thresh)[1,]

    ## 95% Credible Set file
    ### Same as "data.snp" file (as far I can tell), just different format,
    ### and only give the configuration(s) that meet the 95% probability threshold
    # top_CS <- data.table::fread(file.path(locus_dir,"FINEMAP/data.cred"), nThread=nThread)
    # top_CS.snps <- as.character(data.frame(top_CS)[1,][,grep("cred_set_*",colnames(top_CS))])
    # top_CS.probs <- as.numeric(data.frame(top_CS)[1,][,grep("prob_set_*",colnames(top_CS))])
  }
  # CS <- strsplit(top_config$config, ",")[[1]]
  # Import snp-level results
  snp_level <- data.table::fread(file.path(locus_dir,"FINEMAP/data.snp"), nThread = nThread)
  snp_level <- subset(snp_level, prob>credset_thresh & prob_group>credset_thresh) %>%
    dplyr::mutate(CS=1)

  # Merge with original data
  # subset_DT$CS <- ifelse(subset_DT$SNP %in% CS, 1, 0)

  subset_DT <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                             data.table::data.table(subset(snp_level, select=c("rsid","prob","CS")) ),
                                             by.x = "SNP", by.y="rsid")
  subset_DT <- subset_DT %>%
    dplyr::rename(PP=prob) %>%
    dplyr::arrange(dplyr::desc(CS))
  return(subset_DT)
}




#' Retrieve location of \code{FINEMAP} executable
#' @family FINEMAP
#' @keywords internal
#' @source
#' \url{http://www.christianbenner.com}
#' @examples
#' FINEMAP_path <- FINEMAP.find_executable()
FINEMAP.find_executable <- function(FINEMAP_path=NULL,
                                    OS=NULL){
  if(is.null(OS)){OS <- get_os()}
  if(OS=="osx"){exec <- "finemap_v1.3_MacOSX"}else{exec <- "finemap_v1.3.1_x86_64"}
  if(is.null(FINEMAP_path)){
    FINEMAP_path <- system.file("tools",file.path("FINEMAP",exec), package="echolocatoR")
    # FINEMAP_path <- file.path(find.package('echolocatoR'),"exec/FINEMAP/finemap_v1.3_MacOSX")
  }
  return(FINEMAP_path)
}



#' Fine-map locus with \code{FINEMAP}
#'
#' The stepwise conditional search starts with a causal configuration containing the
#' SNP with the lowest P-value alone and then iteratively adds to the causal configuration
#' the SNP given the highest posterior model probability until no further SNP yields
#' a higher posterior model probability.
#'
#' @source
#' \url{http://www.christianbenner.com}
#' @family FINEMAP
#' @keywords internal
#' @examples
#' data("locus_dir"); data("BST1"); data("BST1_LD_matrix");
#' finemap_DT <- BST1
#' dir.create(file.path(locus_dir,"FINEMAP"), showWarnings = FALSE, recursive = TRUE)
#' out <- subset_common_snps(LD_matrix=LD_matrix, finemap_DT=finemap_DT)
#' LD_matrix <- out$LD
#' finemap_DT <- out$DT
#' finemap_DT <- FINEMAP(subset_DT=finemap_DT, locus_dir=locus_dir, LD_matrix=LD_matrix)
FINEMAP <- function(subset_DT,
                    locus_dir,
                    LD_matrix,
                    FINEMAP_path=NULL,
                    n_samples=NULL,
                    n_causal=5,# Max number of allowed causal SNPs
                    model="cond", # cond (stepwise conditional search) vs. sss (stochastic shotgun search)
                    remove_tmps=T,
                    credset_thresh=.95,
                    finemap_version="1.3.1",
                    server=F){
  # n_causal=5; model="cond"; credset_thresh=.95;
  n_samples <- max(subset_DT$N)
  # Setup files
  master_path <- FINEMAP.construct_master(locus_dir = locus_dir,
                                          n_samples = n_samples)
  dat_paths <- FINEMAP.construct_data(locus_dir = locus_dir,
                                      subset_DT = subset_DT,
                                      LD_matrix = LD_matrix)
  # Command line
  ## Example:
  ## cmd <- paste(FINEMAP_path," --sss --in-files",file.path(dirname(FINEMAP_path),"example","master"), "--dataset 1 --n-causal-snps 5")
  # if(startsWith(getwd(),"/sc/")){server <- T}
  # if(server){
  #   FINEMAP_path <- paste0("ml finemap/",finemap_version," && finemap")
  # } else {
    FINEMAP_path <- FINEMAP.find_executable()
    finemap_version <- "1.3"
  # }

  cmd <- paste("cd",locus_dir,"&&",
               FINEMAP_path,
               paste0("--",model),
               "--in-files", master_path,
               "--log",
               # Option to set the maximum number of allowed causal SNPs
               # (Default is 5)
               "--n-causal-snps",n_causal)
  printer(cmd)
  system(cmd)
  # Process results
  finemap_dat <- FINEMAP.process_results(locus_dir = locus_dir,
                                        subset_DT = subset_DT,
                                        credset_thresh = credset_thresh,
                                        finemap_version = finemap_version)
  # Remove tmp files
  if(remove_tmps){
    printer("+FINEMAP: Removing tmp files...")
    tmp_files <- file.path(locus_dir,"FINEMAP",
                           c("data.snp",
                             "data.config",

                             "data.ld",
                             "data.log_cond",
                             "data.log_sss",
                             "data.z",
                             "master")
                           )
    tmp_bool <- suppressWarnings(file.remove(tmp_files))
    tmp_bool <- suppressWarnings(file.remove(file.path(locus_dir,"FINEMAP")))
  }
  return(finemap_dat)
}





