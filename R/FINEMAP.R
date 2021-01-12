# ***************** #
#----  FINEMAP  ----#
# ***************** #


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
#' out <- subset_common_snps(BST1_LD_matrix, finemap_DT)
#' LD_matrix <- out$LD
#' subset_DT <- out$DT
#' subset_DT $N<- subset_DT$N_cases+subset_DT$N_controls
#' finemap_DT <- FINEMAP(subset_DT=subset_DT, locus_dir=locus_dir, LD_matrix=LD_matrix, n_samples=)
FINEMAP <- function(subset_DT,
                    locus_dir,
                    LD_matrix,
                    FINEMAP_path=NULL,
                    n_samples=NULL,
                    n_causal=5,# Max number of allowed causal SNPs
                    model="cond", # cond (stepwise conditional search) vs. sss (stochastic shotgun search)
                    remove_tmps=T,
                    credset_thresh=.95,
                    finemap_version="1.4",
                    server=F,
                    verbose=T){
  # n_causal=5; model="cond"; credset_thresh=.95; verbose=T; finemap_version="1.4"; n_samples=NULL;

  n_samples <- if(is.null(n_samples)) max(subset_DT$N) else n_samples
  dir.create(locus_dir, showWarnings = F, recursive = T)
  # Setup files
  master_path <- FINEMAP.construct_master(locus_dir = locus_dir,
                                          n_samples = n_samples)
  dat_paths <- FINEMAP.construct_data(locus_dir = locus_dir,
                                      subset_DT = subset_DT,
                                      LD_matrix = LD_matrix)
  # Command line
  ## Example:
  ## cmd <- paste(FINEMAP_path," --sss --in-files",file.path(dirname(FINEMAP_path),"example","master"), "--dataset 1 --n-causal-snps 5")
  FINEMAP_path <- FINEMAP.find_executable(version = finemap_version,
                                          verbose = verbose)

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
    printer("+ FINEMAP:: Removing tmp files...")
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



#' Prepare input files for \code{FINEMAP}
#'
#' Creates and saves 1) the summary stats file, and 2) the LD matrix.
#' "Columns beta and se are required for fine-mapping.
#' Column maf is needed to output posterior effect size estimates on the
#' allelic scale. All other columns are not required for computations and
#' can be specified arbitrarily."
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
                                   nThread=4,
                                   verbose=T){
  ####### data.z #######
  if(!"A1" %in% colnames(subset_DT)) {subset_DT$A1 <- "A"; printer("+ FINEMAP:: Optional A1 col missing. Replacing with all 'A's.")};
  if(!"A2" %in% colnames(subset_DT)) {subset_DT$A2 <- "T";  printer("+ FINEMAP:: Optional A2 col missing. Replacing with all 'T's.")};
  if(!"MAF" %in% colnames(subset_DT)) {subset_DT$MAF <- .1; printer(" + FINEMAP:: Optional MAF col missing. Replacing with all '.1's")};
  printer("++ Preparing data.z file for FINEMAP",v=verbose)
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
  printer("++ FINEMAP:: Formatting LD Matrix",v=verbose)
  ## The order of the SNPs in the dataset.ld must correspond to the order of variants in dataset.z.
  # load(file.path(locus_dir,"plink","LD_matrix.RData"))

  # Filter
  data.z <- subset(data.z, rsid %in% rownames(LD_matrix))
  ## This filters AND sorts LD_matrix by the order of rsids in data.z
  LD_filt <- LD_matrix[data.z$rsid, data.z$rsid]

  # Write files
  ## MUST be space-delimited
  printer("++ FINEMAP:: Writing z and ld files...",v=verbose)
  if( dim(data.z)[1]==dim(LD_filt)[1] ){
    # data.z
    data.z_path <- file.path(locus_dir,"FINEMAP","data.z")
    data.table::fwrite(data.z, data.z_path, sep = " ",
                       nThread = 1)
    # Sys.chmod(data.z_path, "777", use_umask = FALSE)
    # data.ld
    data.ld_path <- file.path(locus_dir,"FINEMAP","data.ld")
    data.table::fwrite(data.table:::as.data.table.matrix(LD_filt),
                       data.ld_path, sep=" ", quote = F, col.names = F,
                       nThread = 1)
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
                                     file.k=NA,
                                     verbose=T){ # [optional input]){
  printer("++ FINEMAP:: Constructing master file.",v=verbose)
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
                                    finemap_version="1.4",
                                    results_file=".cred",
                                    nThread=1,
                                    sort_by_CS=T,
                                    verbose=T){
  #### Notes on FINEMAP output files ####
  ##
  ## .snp and .cred are often similiar, but not identical.
  ## Reccomendation: use the .cred file that shows the largest posterior probability for the number of causal variants in line 1 of the file.
  ## and extract credible sets from that file.

  ## Example locus:
  # locus_dir="~/Desktop/Fine_Mapping/Data/GWAS/Marioni_2018/ACE"
  # subset_DT <- data.table::fread(file.path(locus_dir, "Multi-finemap/ACE.Marioni_2018.1KGphase3.multi_finemap.csv.gz"))

  FINEMAP.import_data.snp <- function(locus_dir,
                                      credset_thresh=.95,
                                      verbose=T){
    # NOTES:
    ## .snp files: Posterior probabilities in this file are the marginal posterior probability
    ## that a given variant is causal.
    printer(" + FINEMAP:: Importing marginal probabilties...", v=verbose)
    data.snp <- data.table::fread(file.path(locus_dir,"FINEMAP/data.snp"), nThread = 1)
    data.snp <- subset(data.snp, prob>credset_thresh & prob_group>credset_thresh) %>%
      dplyr::mutate(CS=1)
    return(data.snp)
  }

  FINEMAP.import_data.cred <- function(locus_dir,
                                       verbose=T){
    # NOTES:
    ## .cred files: Conditional posterior probabilities that a given variant is causal
    ## conditional on the other causal variants in the region.
    printer(" + FINEMAP:: Importing conditional probabilties...", v=verbose)
    cred_path <- file.path(locus_dir,"FINEMAP/data.cred")
    data.cred <- data.table::fread(cred_path,
                                   na.strings = c("<NA>","NA"),
                                   nThread = 1)
    cred.cols <- grep("cred*", colnames(data.cred), value = T)
    prob.cols <- grep("prob*", colnames(data.cred), value = T)
    # Restructure data to SNP-wise table format
    CS <- lapply(1:nrow(data.cred), function(i){
      rsids <- subset(data.cred, select=cred.cols)[i,]
      PP_vals <- subset(data.cred, select=prob.cols)[i,]
      cred_sets <- data.table::data.table(SNP=unname( t(rsids)[,1] ),
                                          PP=unname(t(PP_vals)[,1]),
                                          CS=i)
      return(cred_sets)
    }) %>% data.table::rbindlist() %>%
      subset(!is.na(SNP))
    return(CS)
  }

  FINEMAP.import_data.config <- function(locus_dir,
                                         credset_thresh=.95,
                                         pvalue_thresh=.05,
                                         top_config_only=T,
                                         verbose=T){
    # NOTES
    ## .config files: Gives all model results for all the configurations tested
    ## (regardless of whether they're over the 95% probability threshold)
    printer(" + FINEMAP:: Importing config file...", v=verbose)
    config_path <- file.path(locus_dir,"FINEMAP/data.config")
    data.config <- data.table::fread(config_path, nThread=1)

    if(top_config_only){
      data.config <- data.config[1,]
    }
    # Gaurd against future renaming of columns
    if(!is.null(credset_thresh) & "prob" %in% colnames(data.config)){
      data.config <- subset(data.config, prob>=credset_thresh)
    }
    # Not all FINEMAP versions seem to have this "pvalue" column?
    if(!is.null(pvalue_thresh) & "pvalue" %in% colnames(data.config)){
      data.config <- subset(data.config, pvalue<pvalue_thresh)
    }

    # Restructre config file
    ## Use the probability of the configuration itself as the snp-wise probabilties
    data.config_format <- data.frame(SNP=strsplit(data.config$config, ",")[[1]],
                                     PP=data.config$prob,
                                     CS=1)
    return(data.config_format)
  }


  #### Process FINEMAP results ####
  if(results_file==".cred"){
    dat <- FINEMAP.import_data.cred(locus_dir = locus_dir,
                                    verbose = verbose)
    # Merge with original dataframe
    subset_DT <- data.table::merge.data.table(data.table::data.table(subset_DT),
                                              data.table::data.table(dat),
                                              by="SNP",
                                              all.x = T)
  }
  if (results_file==".snp"){
    dat <- FINEMAP.import_data.snp(locus_dir = locus_dir,
                                   verbose = verbose) %>%
      dplyr::rename(PP=prob)
    # Merge with original dataframe
    subset_DT <- data.table::merge.data.table(data.table::data.table(subset_DT),
                                              data.table::data.table(subset(dat, select=c("rsid","prob","CS")) ),
                                              by.x = "SNP",
                                              by.y="rsid",
                                              all.x = T)
  }
  if (results_file==".config"){
    dat <- FINEMAP.import_data.config(locus_dir = locus_dir,
                                      verbose = verbose)
    subset_DT <- data.table::merge.data.table(data.table::data.table(subset_DT),
                                              data.table::data.table(dat),
                                              by="SNP",
                                              all.x = T)
  }
  # Sort so that CS SNPs are at the top
  if(sort_by_CS){
    subset_DT <- subset_DT %>% dplyr::arrange(dplyr::desc(PP))
  }
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
                                    OS=NULL,
                                    version="1.4",
                                    verbose=T){
  if(is.null(OS)){OS <- get_os()}

  if(version=="1.4"){
    printer("+ Using FINEMAP v1.4",v=verbose)
    if(OS=="osx"){
      exec <- "finemap_v1.4_MacOSX"
    } else{
      exec <- "finemap_v1.4_x86_64"
    }
  } else {
    printer("+ Using FINEMAP v1.3",v=verbose)
    if(OS=="osx"){
      exec <- "finemap_v1.3_MacOSX"
    } else{
      exec <- "finemap_v1.3.1_x86_64"
    }
  }

  if(is.null(FINEMAP_path)){
    FINEMAP_path <- system.file("tools",file.path("FINEMAP",exec), package="echolocatoR")
    # FINEMAP_path <- file.path(find.package('echolocatoR'),"exec/FINEMAP/finemap_v1.3_MacOSX")
  }
  return(FINEMAP_path)
}






