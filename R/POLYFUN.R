# $$$$$$$$$$$$$$$$ $$$$$$$ $$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$ PolyFun $$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$ $$$$$$$ $$$$$$$$$$$$$$$$
# https://github.com/omerwe/polyfun

##-------------------------------------------------------------
# There are three ways to run PolyFun:

# 1. Using precomputed prior causal probabilities of 19 million imputed UK Biobank SNPs with MAF>0.1%, based on a meta-analysis of 15 UK Biobank traits. This is the simplest approach, but it may not include all your SNPs of interest (especially when analyzing non-European populations) and the prior causal probabilities may not be optimal for some traits.

# 2. Computing prior causal probabilities via an L2-regularized extension of stratified LD-score regression (S-LDSC). This is a relatively simple approach, but the prior causal probabilities may not be robust to modeling misspecification.

# 3. Computing prior causal probabilities non-parametrically. This is the most robust approach, but it is computationally intensive and requires access to individual-level genotypic data from a large reference panel (optimally >10,000 population-matched individuals).
##-------------------------------------------------------------


# GitHub Notes:
## How to pull changes from original repo into the forked repo
# https://digitaldrummerj.me/git-syncing-fork-with-original-repo/



#' Read parquet file as data.frame
#'
#' @param method Specify whether to use \code{\link{SparkR}} (R)
#' or \code{pandas} (Python + \code{reticulate}).
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' parquet_path <- system.file("tools/polyfun/example_data/weights.10.l2.ldscore.parquet", package = "echolocatoR")
#' # Using python (pandas) - default
#' dat <- POLYFUN.read_parquet(parquet_path=parquet_path, method="pandas")
#' # Using R (SparkR)
#' dat <- POLYFUN.read_parquet(parquet_path=parquet_path, method="sparkR")
#' }
POLYFUN.read_parquet <- function(parquet_path,
                                 conda_env="echoR",
                                 method="pandas"){
  if(method=="SparkR" & "SparkR" %in% rownames(utils::installed.packages())){
    printer("+ Importing parquet file with `SparkR (R)`")
    SparkR::sparkR.session()
    parquor <- SparkR::read.parquet(parquet_path)
    parquor <- SparkR::as.data.frame(parquor) %>%
      data.table::data.table()
  } else {
    printer("+ Importing parquet file with `pandas (Python)`")
    python <- CONDA.find_python_path(conda_env = conda_env)
    Sys.setenv(RETICULATE_PYTHON=python)
    reticulate::use_python(python = python)
    CONDA.activate_env(conda_env = conda_env)
    # pandas <- CONDA.find_package(package = "pandas", conda_env = conda_env)
    pd <- reticulate::import("pandas")
    parquor <- pd$read_parquet(parquet_path)
    system(paste(python))
  }
  return(parquor)
}




#' Display PolyFun help
#' @keywords internal
#' @family polyfun
#' @examples
#' POLYFUN.help()
POLYFUN.help <- function(polyfun=NULL,
                         conda_env="echoR"){
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)
  python <- CONDA.find_python_path(conda_env = conda_env)
  cmd <- paste(python,
               file.path(polyfun,"polyfun.py"),
               "--help")
  system(cmd)
}




#' Folder where PolyFun submodule is stored
#' @keywords internal
#' @family polyfun
#' @examples
#' polyfun_path <- POLYFUN.find_polyfun_folder()
POLYFUN.find_polyfun_folder <- function(polyfun_path=NULL){
  if(is.null(polyfun_path)){
    polyfun_path <- system.file("tools/polyfun","",package = "echolocatoR")
  }
  return(polyfun_path)
}





# %%%%%%%%%%%%%%%% PolyFun approach 1 %%%%%%%%%%%%%%%%
## Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Biobank traits



#' Create output dir and import SNP data.frame
#' @keywords internal
#' @family polyfun
#' @examples
#' data("BST1"); data("locus_dir");
#' finemap_DT <- POLYFUN.initialize(locus_dir=locus_dir, finemap_dat=BST1)
POLYFUN.initialize <- function(locus_dir,
                               finemap_dat=NULL,
                               nThread=4){
  dataset <- basename(dirname(locus_dir))
  locus <- basename(locus_dir)
  # Create path
  PF.output.path <- file.path(locus_dir, "PolyFun")
  dir.create(PF.output.path, showWarnings = F, recursive = T)
  # Import SNPs
  if(is.null(finemap_dat)){
    if(is.null(locus)){
      printer("POLYFUN:: Importing summary stats from disk: genome-wide")
      finemap_dat <- data.table::fread(Directory_info(dataset, "fullSS.local"), nThread = 4)
    }
    printer("POLYFUN:: Importing summary stats from disk:",locus)
    finemap_dat <- data.table::fread(file.path(dirname(Directory_info(dataset, "fullSS.local")),locus,
                                              paste(locus,dataset,"subset.tsv.gz",sep="_")), nThread = nThread)
  }
  return(finemap_dat)
}




#' Prepare SNP input for PolyFun
#'
#' PolyFun requires a space-delimited (gzipped or not) file with these columns:
#' \itemize{
#' \item{CHR}
#' \item{BP}
#' \item{A1}
#' \item{A2}
#' }
#' @keywords internal
#' @family polyfun
#' @examples
#' data("BST1"); data("locus_dir");
#' finemap_dat <- BST1
#' PF.output.path <- file.path(locus_dir, "PolyFun")
#' POLYFUN.prepare_snp_input(PF.output.path=PF.output.path, locus_dir=locus_dir, finemap_dat=finemap_dat)
POLYFUN.prepare_snp_input <- function(PF.output.path,
                                      locus_dir,
                                      finemap_dat=NULL,
                                      nThread=1){
  printer("PolyFun:: Preparing SNP input file...")
  if(!"A1" %in% colnames(finemap_dat)) A1 <- NULL;
  if(!"A2" %in% colnames(finemap_dat)) A2 <- NULL;
  PF.dat <- dplyr::select(finemap_dat,
                          SNP,
                          CHR,
                          BP=POS,
                          A1=A1,
                          A2=A2)
  printer("+ PolyFun::",nrow(PF.dat),"SNPs identified.")
  snp.path <- file.path(PF.output.path,"snps_to_finemap.txt.gz")
  printer("+ PolyFun:: Writing SNP file ==>",snp.path)
  dir.create(dirname(snp.path), recursive = T, showWarnings = F)
  data.table::fwrite(PF.dat, file = snp.path,
                     nThread = nThread, sep = " ")
  # Check if it's actually gzipped
  ## (in older versions of data.table, files with the .gz extension are automatically gzipped)
  if(!R.utils::isGzipped(snp.path)){
    R.utils::gzip(filename=snp.path, destname=snp.path, overwrite=T)
  }
  return(snp.path)
}








#' Extract pre-computed prior
#'
#' Extract SNP-wise prior probabilities pre-computed from many UK Biobank traits.
#'
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#'
#' data("BST1"); data("locus_dir");
#' priors <- POLYFUN.get_precomputed_priors(locus_dir=locus_dir, finemap_dat=BST1)
#' }
POLYFUN.get_precomputed_priors <- function(polyfun=NULL,
                                           locus_dir,
                                           finemap_dat=NULL,
                                           force_new_priors=T,
                                           remove_tmps=F,
                                           nThread=4,
                                           conda_env="echoR"){
  python <- CONDA.find_python_path(conda_env = conda_env)
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)
  dataset <- basename(dirname(locus_dir))
  locus <- basename(locus_dir)
  PF.output.path <- file.path(locus_dir, "PolyFun")
  snp_w_priors.file <- file.path(PF.output.path,"snps_with_priors.snpvar.tsv.gz")

  if((file.exists(snp_w_priors.file)) & force_new_priors==F){
    print("++ Importing pre-existing priors.")
    priors <- data.table::fread(snp_w_priors.file, nThread = nThread) %>%
      dplyr::rename(SNP=SNP_x) %>% dplyr::select(-SNP_y)
    return(priors)
  } else {
    finemap_dat <- POLYFUN.initialize(finemap_dat=finemap_dat,
                                     locus_dir=locus_dir)
    chrom <- finemap_dat$CHR[1]
    # [1]. Prepare input
    snp.path <- POLYFUN.prepare_snp_input(PF.output.path=PF.output.path,
                                          locus_dir=locus_dir,
                                          finemap_dat=finemap_dat)
    # [2.] Retrieve priors
    # test <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/LRRK2/PolyFun/snps_to_finemap.txt.gz")
    try({
      cmd <- paste(python,
                   file.path(polyfun,"extract_snpvar.py"),
                   "--snps",snp.path,
                   "--out",snp_w_priors.file)
      print(cmd)
      system(cmd)
      printer("++ Remove tmp file.")
    })

    miss.file <- paste0(snp_w_priors.file,".miss.gz")
    if(file.exists(miss.file)){
      printer("+ PolyFun:: Rerunning after removing missing SNPs.")
      miss.snps <- data.table::fread(miss.file)
      filt_DT <- subset(finemap_dat, !(SNP %in% miss.snps$SNP) )
      if(remove_tmps){file.remove(miss.file)}
      priors <- POLYFUN.get_precomputed_priors(locus_dir=locus_dir,
                                               finemap_dat=filt_DT,
                                               force_new_priors=force_new_priors,
                                               conda_env = conda_env)
    }
    # Import results
    priors <- data.table::fread(snp_w_priors.file,
                                nThread = nThread) %>%
      dplyr::rename(SNP=SNP_x) %>% dplyr::select(-SNP_y)
    return(priors)
  }
  if(remove_tmps){ file.remove(snp.path) }
}




# %%%%%%%%%%%%%%%% PolyFun approaches 2 & 3 %%%%%%%%%%%%%%%%
## 2. Computing prior causal probabilities via an L2-regularized extension of S-LDSC
## 3. Computing prior causal probabilities non-parametrically


#' Munge summary stats
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' data("genome_wide_dir");
#' fullSS_path <- example_fullSS(fullSS_path="~/Desktop/Nalls23andMe_2019.fullSS_subset.tsv")
#' munged_path <- POLYFUN.munge_summ_stats(fullSS_path=fullSS_path, locus_dir=genome_wide_dir, force_new_munge=T)
#' }
POLYFUN.munge_summ_stats <- function(polyfun=NULL,
                                     fullSS_path,
                                     locus_dir,
                                     sample_size=NULL,
                                     min_INFO=0,
                                     min_MAF=0.001,
                                     force_new_munge=F,
                                     conda_env="echoR"){
  python <- CONDA.find_python_path(conda_env = conda_env)
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)
  PF.output.path <- file.path(locus_dir, "PolyFun")
  dir.create(PF.output.path, showWarnings = F, recursive = T)
  munged_path <- file.path(PF.output.path,
                           paste0(gsub("\\.gz|\\.txt|\\.tsv|\\.csv","",basename(fullSS_path)),".munged.parquet"))
  # PolyFun requires space-delimited file with the following columns (munging can recognize several variations of these names):
  ## SNP CHR BP ....and....
  ## either a p-value, an effect size estimate and its standard error, a Z-score or a p-value

  header <- get_header(large_file = fullSS_path)
  if("N_cases" %in% header & "N_controls" %in% header){
    warning("Cannot both specify sample_size (--n) and have an N_cases/N_controls column in the sumstats file. Using N_cases/N_controls columns instead.")
    sample_size_arg <- NULL
  } else {sample_size_arg <- paste("--n",sample_size)}

  if(!file.exists(munged_path) | force_new_munge){
    printer("+ PolyFun:: Initiating data munging pipeline...")
    cmd <- paste(python,
                 file.path(polyfun,"munge_polyfun_sumstats.py"),
                 "--sumstats",fullSS_path,#Directory_info(dataset_name = dataset, ifelse(server,"fullSS",fullSS.loc)),
                 sample_size_arg, # Study sample size
                 "--out",munged_path,
                 "--min-info",min_INFO,
                 "--min-maf",min_MAF
                 )
    print(cmd)
    system(cmd)
  } else {printer("+ PolyFun:: Existing munged summary stats files detected.")}
  return(munged_path)
}




#' Find and import LD score files
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' output_prefix <- file.path(system.file("tools/polyfun/gold","",package = "echolocatoR"),"testrun.22")
#' ldscore <- POLYFUN.gather_ldscores(output_prefix=output_prefix)
#' }
POLYFUN.gather_ldscores <- function(output_prefix){
  ldscore.files <-  list.files(dirname(output_prefix),
                               pattern = ".l2.ldscore.parquet", full.names = T)
  if(length(ldscore.files)>1){printer("POLYFUN:: >1 ldscore file detected. Only using the first:",ldscore.files[1])}
  parquor <- POLYFUN.read_parquet(ldscore.files[1])
  return(parquor)
}




#' Find and import PolyFun annotation files
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' data("BST1")
#' finemap_DT <- BST1
#' subset_snps <- finemap_DT$SNP
#' annot_DT <- POLYFUN.gather_annotations(chromosomes=finemap_DT$CHR[1], subset_snps=subset_snps, polyfun_annots="/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB")
#' }
POLYFUN.gather_annotations <- function(chromosomes=c(1:22),
                                       subset_snps=NULL,
                                       polyfun_annots){
  annot_DT <- lapply(chromosomes, function(chrom){
    printer("+ POLYFUN:: Chromosome",chrom,"...")
    parquet.file <- list.files(path = polyfun_annots,
                               # e.g. " /pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/baselineLF2.2.UKB.5.annot.parquet"
                               pattern = paste0("*\\.",chrom,"\\.annot\\.parquet"),
                               full.names = T)
    annot_df <- POLYFUN.read_parquet(parquet_path = parquet.file)
    annot_df$BP <- as.integer(annot_df$BP)
    if(!is.null(subset_snps)){
      annot_df <- subset(annot_df, SNP %in% subset_snps)
    }
    return(data.table::data.table(annot_df))
  }) %>% data.table::rbindlist()
  return(annot_DT)
}




#' Download 1000 Genomes reference files
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' ref.prefix <- POLYFUN.download_ref_files(force_overwrite=T)
#' }
POLYFUN.download_ref_files <- function(alkes_url="https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_plinkfiles.tgz",
                                       # output_path="/sc/arion/projects/pd-omics/data/1000_Genomes/Phase1",
                                       results_dir="./results",
                                       force_overwrite=F,
                                       download_method="wget"
                                       ){
  output_path <- file.path(results_dir,"resources/1000Genomes_Phase1")
  file_name <- basename(alkes_url)
  dir.create(output_path, showWarnings = F, recursive = T)
  # cmd <- paste("wget",
  #              alkes_url,
  #              "--no-parent", # Makes everything waayyyyyy faster
  #              "& mv",file_name, output_path,
  #              "& tar zxvf",output_path,"--strip 1")
  # paste(cmd)
  # system(cmd)
  printer("+ POLYFUN:: Downloading reference files...")
  out_file <- downloader(input_url = alkes_url,
                         output_path = output_path,
                         force_overwrite = force_overwrite,
                         download_method = download_method,
                         background = F)
  printer("+ POLYFUN:: Unzipping reference files...")
  system(paste("tar zxvf",file.path(output_path, file_name),
               "--strip 1",
               ifelse(force_overwrite,"","-k"),
               "-C",output_path))
  # List reference files
  ref.prefix <- list.files(output_path, pattern = "*.bim", full.names = T)[1]
  ref.prefix <- gsub(".?.bim","", ref.prefix)
  return(ref.prefix)
}





#' Get ref data path prefix
#'
POLYFUN.get_ref_prefix <- function(ref_prefix=NULL,
                                   locus_dir){
  "./resources/1000_Genomes/Phase1/1000G.mac5eur."


}

#' Recompute SNP-wise priors from summary stats
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' data("locus_dir");
#' fullSS_path <- example_fullSS(fullSS_path="~/Desktop/Nalls23andMe_2019.fullSS_subset.tsv")
#' munged_path <- "results/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/nallsEtAl2019_allSamples_allVariants.mod.munged.parquet"
#' LDSC.files <- POLYFUN.compute_priors(locus_dir=locus_dir, fullSS_path=fullSS_path, conda_env="echoR")
#' }
POLYFUN.compute_priors <- function(polyfun=NULL,
                                    locus_dir,
                                    fullSS_path,
                                    sample_size = NULL,
                                    min_INFO = 0,
                                    min_MAF = 0.001,
                                    annotations_path=NULL,
                                    weights_path=NULL,
                                    prefix="PD_GWAS",
                                    chrom="all",
                                    compute_ldscores=F,
                                    allow_missing_SNPs=T,
                                    ref_prefix=NULL,
                                    remove_tmps=T,
                                    conda_env = "echoR"){
  # polyfun="./echolocatoR/tools/polyfun"; parametric=T;  weights.path=file.path(polyfun,"example_data/weights."); annotations.path=file.path(polyfun,"example_data/annotations."); munged.path= "./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/sumstats_munged.parquet"; parametric=T; dataset="Nalls23andMe_2019"; prefix="PD_GWAS"; compute_ldscores=F; allow_missing_SNPs=T; chrom="all"; finemap_dat=NULL; locus="LRRK2"; server=F; ref.prefix="/sc/arion/projects/pd-omics/data/1000_Genomes/Phase1/1000G.mac5eur.";
  # CONDA.activate_env(conda_env = conda_env)
  # PATH_cmd <- "source ~/.bash_profile &&"
  python <- CONDA.find_python_path(conda_env = conda_env)
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)



  if(is.null(annotations_path)){annotations_path <- file.path(system.file("tools/polyfun/example_data",package="echolocatoR"),"annotations.")}
  if(is.null(weights_path)){weights_path <- file.path(system.file("tools/polyfun/example_data",package="echolocatoR"),"weights.")}
  # 0. Create paths
  PF.output.path <- file.path(locus_dir, "PolyFun")
  dir.create(PF.output.path, showWarnings = F, recursive = T)
  out.path <- file.path(PF.output.path,"output")
  output_prefix <- file.path(out.path, prefix, prefix)
  dir.create(out.path, showWarnings = F, recursive = T)


  # 1. Munge summary stats

  printer("PolyFun:: [1]  Create a munged summary statistics file in a PolyFun-friendly parquet format.")
  munged.path <- POLYFUN.munge_summ_stats(polyfun=polyfun,
                                           fullSS_path = fullSS_path,
                                           locus_dir=locus_dir,
                                           sample_size=sample_size,
                                           min_INFO = min_INFO,
                                           min_MAF = min_MAF,
                                           force_new_munge = F,
                                           conda_env = conda_env)

  # 2.
  ## If compute_ldscores == F:
  # This will create 2 output files for each chromosome: output/testrun.<CHR>.snpvar_ridge.gz and output/testrun.<CHR>.snpvar_ridge_constrained.gz. The first contains estimated per-SNP heritabilities for all SNPs (which can be used for downstream analysis with PolyFun; see below), and the second contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping.
  # library(reticulate)
  # reticulate::use_virtualenv("echoR")
  # pd <- reticulate::import("pandas")
  # pd$read_csv("./Data/directories_table.csv")
  # reticulate::
  # source_python(file.path(polyfun,"polyfun.py"))

  # NOTE! if you're running without the "--no-partitions" flag,
  ## you need to load R first `ml R`.
  printer("PolyFun:: [2] Run PolyFun with L2-regularized S-LDSC")
  # pf <- reticulate::py_run_file(polyfun)
  # reticulate::py_call()

  cmd2 <- paste(python,
                file.path(polyfun,"polyfun.py"),
                "--compute-h2-L2",
               # Approach 2 = Parametric = no partitions = T
               # Approach 3 = Non-parametric = partitions = F
                ifelse(compute_ldscores,"","--no-partitions"),
                "--output-prefix",output_prefix,
                "--sumstats",munged_path,
                "--ref-ld-chr",annotations_path,
                "--w-ld-chr",weights_path,
                ifelse(allow_missing_SNPs,"--allow-missing",""))
  print(cmd2)
  system2(cmd2)

  # Computationally intensive: can parallelize by chromosomes
  if(compute_ldscores){
    # 3. Computationally intensive step
    printer("PolyFun:: [3] Compute LD-scores for each SNP bin")
    cmd3 <- paste(python,
                  file.path(polyfun,"polyfun.py"),
                  "--compute-ldscores",
                  "--output-prefix",output_prefix,
                  "--bfile-chr",ref_prefix,
                  ifelse(chrom=="all","",paste("--chr",chrom)),
                  ifelse(allow_missing_SNPs,"--allow-missing","") )
    print(cmd3)
    system2(cmd3)
    # 4.
    printer("PolyFun:: [4] Re-estimate per-SNP heritabilities via S-LDSC")
    cmd4 <- paste(python,
                  file.path(polyfun,"polyfun.py"),
                  "--compute-h2-bins",
                  "--output-prefix",output_prefix,
                  "--sumstats",munged_path,
                  "--w-ld-chr",weights_path,
                  ifelse(allow_missing_SNPs,"--allow-missing",""))
    print(cmd4)
    system(cmd4)

    printer("PolyFun:: Results directory =",dirname(output_prefix))
    printer("PolyFun:: Results files:")
    printer("          *.snpvar_ridge.gz")
    printer("          *.snpvar_ridge_constrained.gz")
    # The output of the PARTITIONED LDSC has the suffix: .snpvar_constrained.gz (one per chrom)
    LDSC.files <- list.files(out.path,
                             pattern = "*.snpvar_constrained.gz", full.names = T)
    # pd_ldsc <- data.table::fread(PS_LDSC.files[1], nThread = 4)
    # ldscore <- POLYFUN.read_parquet(file.path(out.path,"PD_GWAS.1.l2.ldscore.parquet"))
    # bin.1 <- POLYFUN.read_parquet(file.path(out.path,"PD_GWAS.2.bins.parquet"))
    #rowSums(bin.1[,-c(1:5)]) # each SNP belongs to only 1 bin
  } else { LDSC.files <- list.files(out.path, pattern = "_ridge_constrained.gz", full.names = T) }



  if(remove_tmps){ file.remove(snp.path) }
  return(LDSC.files)
}







#### Original LDSC (extended by PolyFun)



#' Run a modified version of S-LDSC
#'
#' Modifications to S-LDSC include L2-regularization.
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.run_ldsc <- function(polyfun=NULL,
                             output_dir=NULL,
                             munged.path,
                             min_INFO = 0.6,
                             min_MAF = 0.05,
                             annotations.path=file.path(polyfun,"example_data/annotations."),
                             weights.path=file.path(polyfun,"example_data/weights."),
                             prefix="LDSC",
                             chrom="all",
                             compute_ldscores=F,
                             allow_missing_SNPs=T,
                             munged_path="/sc/arion/projects/pd-omics/tools/polyfun/Nalls23andMe_2019.sumstats_munged.parquet",
                             ref.prefix="/sc/arion/projects/pd-omics/data/1000_Genomes/Phase1/1000G.mac5eur.",
                             freq.prefix="/sc/arion/projects/pd-omics/tools/polyfun/1000G_frq/1000G.mac5eur.",
                             conda_env="echoR"){
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)
  python <- CONDA.find_python_path(conda_env = conda_env)
  # if(server){
  #   annotations.path <-  "/sc/arion/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/baselineLF2.2.UKB."
  #   weights.path <-  "/sc/arion/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/weights.UKB."
  # }

  # 0. Create paths
  dir.create(output_dir, showWarnings = F, recursive = T)
  out.path <- file.path(output_dir,"output")
  output_prefix <- file.path(out.path, prefix, prefix)
  dir.create(out.path, showWarnings = F, recursive = T)
  # https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
  cmd <- paste(python,
               file.path(polyfun,"ldsc.py"),
                "--h2",munged_path,
                "--ref-ld-chr",annotations.path,
                "--w-ld-chr",weights.path,
                "--overlap-annot",
                "--frqfile-chr",freq.prefix,
                "--not-M-5-50", # Important! enrichment estimates will be provided with MAF>0.1% SNPs instead of MAF>5% SNPs.
                "--out",output_prefix)
  # help_cmd <- paste("python",file.path(polyfun,"ldsc.py -h"))
  print(cmd)
  system(cmd)
}






# %%%%%%%%%%%%%%%% Run PolyFun+SUSIE %%%%%%%%%%%%%%%%


#' Run PolyFun+SUSIE fine-mapping pipeline
#'
#' Uses echolocatoR wrapper for SUSIE instead of the \code{\link{POLYFUN.finemapper}}
#' function which uses a python script provided with PolyFun.
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN_SUSIE <- function(locus_dir,
                          polyfun=NULL,
                          finemap_dat=NULL,
                          LD_matrix=NULL,
                          polyfun_approach="non-parametric",
                          dataset_type="GWAS",
                          max_causal=5,
                          sample_size=NULL,
                          server=F,
                          PP_threshold=.95,
                          conda_env="echoR"){
  # polyfun="./echolocatoR/tools/polyfun";  locus_dir="./Data/GWAS/Nalls23andMe_2019/_genome_wide"; dataset="Nalls23andMe_2019"; locus="LRRK2"; finemap_dat=NULL; polyfun_priors="parametric"; sample_size=1474097; min_INFO=0; min_MAF=0; server=T; dataset_type="GWAS"; n_causal=10; PP_threshold=.95
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)
  out.path <- file.path(dirname(locus_dir),"_genome_wide/PolyFun/output")
  chrom <- unique(finemap_dat$CHR)
  # printer("++ POLYFUN:: Unique chrom =",paste(chrom,collapse=","))

  # Import priors
  # ~~~~~~~~ Approach 1 ~~~~~~~~
  if (polyfun_approach=="precomputed"){
    priors <- POLYFUN.get_precomputed_priors(locus_dir=locus_dir,
                                             finemap_dat=finemap_dat,
                                             force_new_priors=F,
                                             conda_env=conda_env)
    # precomputed.priors <- priors
  # ~~~~~~~~ Approach 2 ~~~~~~~~
  } else if (polyfun_approach=="parametric"){
    ldsc.files <- list.files(out.path, pattern = "*.snpvar_ridge_constrained.gz", full.names = T) %>%
      grep(pattern = paste0(".",chrom,"."), value = T, fixed=T)
    priors <- .rbind.file.list(ldsc.files)
    # ~~~~~~~~ Approach 3 ~~~~~~~~
  } else if (polyfun_approach=="non-parametric"){
    ldsc.files <- list.files(out.path, pattern = "*.snpvar_constrained.gz", full.names = T) %>%  base::grep(pattern = paste0(".",chrom,"."), value = T, fixed = T)
    priors <- .rbind.file.list(ldsc.files)
  }

  # Ensure formatting is correct (sometimes SNP gets turned into logical?)
  finemap_dat$SNP <- as.character(finemap_dat$SNP)
  priors <- dplyr::select(priors, SNP, POLYFUN.h2=SNPVAR) %>%
    data.table::data.table() %>%
    dplyr::mutate(SNP=as.character(SNP))
  # Prepare data
  merged_dat <- data.table::merge.data.table(x = finemap_dat,
                                             y = priors,
                                             by="SNP")
  sub.out <- subset_common_snps(LD_matrix = LD_matrix,
                                finemap_dat = merged_dat)
  LD_matrix <- sub.out$LD
  new_DT <- sub.out$DT
  # Run SUSIE
  subset_DT <- SUSIE(subset_DT=new_DT,
                     LD_matrix=LD_matrix,
                     dataset_type=dataset_type,
                     max_causal=max_causal,
                     sample_size=sample_size,
                     PP_threshold=PP_threshold,

                     prior_weights=new_DT$POLYFUN.h2,
                     rescale_priors = T)
  # subset_DT.precomputed <- subset_DT
  # subset_DT.computed <- subset_DT
  # Check for differences between pre-computed and re-computed heritabilities
  # library(patchwork)
  # ## GWAS
  # ggplot() +
  #   geom_point(data=subset_DT, aes(x=POS, y=-log10(P), fill="PD GWAS", color=-log10(P))) +
  #   #3 Priors
  #   ggplot() +
  #   geom_point(data=precomputed.priors, aes(x=BP, y=SNPVAR, color="Pre-computed")) +
  #   geom_point(data=subset(priors, BP>=min(precomputed.priors$BP) & BP<=max(precomputed.priors$BP)),
  #              aes(x=BP, y=SNPVAR, color="Computed")) +
  #   ## fine-mapping with computed priors
  #   ggplot() +
  #   geom_point(data=subset_DT.computed, aes(x=POS, y=PP, color="POLYFUN+SUSIE \n Computed Priors")) +
  #   patchwork::plot_layout(ncol = 1) +
  #   ## fine-mapping with pre-computed priors
  #   ggplot() +
  #   geom_point(data=subset_DT.precomputed, aes(x=POS, y=PP, fill="POLYFUN+SUSIE \n Pre-computed Priors"), color="turquoise")
  return(subset_DT)
}





#' Run \emph{PolyFun+SUSIE} fine-mapping pipeline
#'
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.finemapper <- function(polyfun=NULL,
                               finemap_dat=NULL,
                               npz_gz_LD=NULL,
                               locus=NULL,
                               sample_size=NULL,
                               locus_dir,#="Data/GWAS/Nalls23andMe_2019/LRRK2",
                               n_causal=5,
                               method="susie",
                               h2_path=NULL,
                               conda_env="echoR"){
  # base_url  <- "./echolocatoR/tools/polyfun/LD_temp"
  polyfun <- POLYFUN.find_polyfun_folder(polyfun_path = polyfun)
  python <- CONDA.find_python_path(conda_env = conda_env)

  chrom <- unique(finemap_dat$CHR)
  file.name <- paste0("chr",chrom,"_","40000001_43000001")
  ld_path <- file.path(locus_dir,"LD",file.name)
  out_path <- file.path(locus_dir,"PolyFun",
                        paste0("finemapper_res.",basename(locus_dir),".gz"))
  dir.create(dirname(out_path), showWarnings = F, recursive = T)

  # if(is.null(h2_path)){
  #   url_prefix <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output/PD_GWAS."
  #   h2_path <- paste0(url_prefix,chrom,".snpvar_constrained.gz")
  # }

  if(is.null(sample_size) & (!"N" %in% colnames(finemap_dat))){
    finemap_dat <- get_sample_size(finemap_dat)
    sample_size <- max(finemap_dat$N)
  }
  # munged.path <- "~/Desktop/Nalls23andMe_2019.sumstats_munged.parquet"
  # locus="LRRK2"
  if(is.null(npz_gz_LD)){
    npz_path <- list.files(file.path(locus_dir, "LD"), ".npz$",full.names = T)
    if(length(npz_path)>0){
      npz_path <- npz_path[1]
    } else {
      rds_path <- list.files(file.path(locus_dir, "LD"), ".RDS$",full.names = T)[1]
      npz_path <- LD.rds_to_npz(rds_path = rds_path)
    }
    npz_prefix <- gsub(".npz$","",npz_path)
  }

  cmd <- paste(python,
               file.path(polyfun,"run_finemapper.py"),
                "--ld",npz_prefix,
                "--sumstats",h2_path,
                "--n",sample_size,
                "--chr",chrom,
                "--start",min(finemap_dat$POS),
                "--end",max(finemap_dat$POS),
                "--method",method,
                "--max-num-causal",n_causal,
                # "--threads 2",# use max detected cores if not specified
                "--out",out_path)
  print(cmd)
  system(cmd)

  h2 <- .rbind.file.list("Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output/PD_GWAS.16.snpvar_constrained.gz")
}




#' Plot PolyFun and other fine-mapping results
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.plot <- function(subset_DT,
                         LD_matrix,
                         locus=NULL,
                         subtitle="PolyFun Comparison",
                         conditions=c("SUSIE","POLYFUN_SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")){
  # Quickstart
  # locus="LRRK2"; subtitle="PolyFun Comparison"; conditions=c("SUSIE","POLYFUN_SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")
  # # Get r2

  if(plot_ld){
    lead.snp <- top_n(subset_DT,1,-P)$SNP #subset(subset_DT, leadSNP==T)$SNP
    r2 <- data.table::data.table(SNP=names(LD_matrix[lead.snp,]),
                                 r2=LD_matrix[lead.snp,]^2)
    dat <- data.table:::merge.data.table(subset_DT, r2, by="SNP")
  } else{dat <- dplyr::mutate(subset_DT, r2=1)}
  dat <- dplyr::mutate(dat, Mb=round(POS/1000000,3))



  library(patchwork)
  # GWAS
  gg <- ggplot(dat, aes(x=Mb, y=-log10(P), color=r2)) +
    scale_color_gradient(low="blue",high="red", breaks=c(0,.5,1), limits=c(0,1)) +
    geom_point() +
    labs(y="GWAS -log10(P)") +
    ggrepel::geom_label_repel(data = top_n(dat,n=1,-P),
                              aes(label=SNP),alpha=0.7) +
    geom_point(data=top_n(dat,n=1,-P), size=5, shape=1, color="red") +
    scale_y_continuous(limits = c(0,max(-log10(dat$P))*1.1)) +

    # PolyFun priors
    ggplot(dat, aes(x=Mb, y=POLYFUN.h2, color=POLYFUN.h2)) +
    scale_color_viridis_c(limits=c(0,1), breaks=c(0,.5,1)) +
    geom_point() +
    # ylim(0,1) +

    # PolyFun+SUSIE PP
    ggplot(dat, aes(x=Mb, y=POLYFUN_SUSIE.PP, color=POLYFUN_SUSIE.PP)) +
    geom_point() +
    ggrepel::geom_label_repel(data = subset(dat,POLYFUN_SUSIE.PP>=.5),
                              aes(label=SNP),alpha=0.7, color='green') +
    geom_point(data=subset(dat, POLYFUN_SUSIE.CS>0),
               #subset(dat, PolyFun_SUSIE.PP>=.95),
               size=5, shape=1, color="green") +
    scale_y_continuous(breaks = c(0,.5,1), limits = c(0,1.1)) +
    scale_color_continuous(breaks=c(0,.5,1), limits=c(0,1)) +

    # SUSIE PP
    ggplot(dat, aes(x=Mb, y=SUSIE.PP, color=SUSIE.PP)) +
    geom_point() +
    ggrepel::geom_label_repel(data = subset(dat,SUSIE.CS>0),
                              aes(label=SNP),alpha=0.8, color='green') +
    geom_point(data=subset(dat,SUSIE.PP>=.95),
               size=5, shape=1, color="green") +
    scale_y_continuous(breaks = c(0,.5,1), limits = c(0,1.1)) +
    scale_color_continuous(breaks=c(0,.5,1), limits=c(0,1)) +
    # Overall layers
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(title = locus,
                               subtitle = paste(nrow(dat),"SNPs"),#"PolyFun Comparison",
                               theme =  theme(plot.title = element_text(hjust = 0.5),
                                              plot.subtitle = element_text(hjust = 0.5)))
  print(gg)
  ggsave(file.path(locus_dir,'PolyFun',"PolyFun.plot.png"), plot = gg,
         dpi = 400, height = 10, width = 7)
}





## ------------------------------------------ ##
########## POLYFUN H2 ENRICHMENT ##############
## ------------------------------------------ ##




#' Run and plot heritability enrichment tests
#'
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.ldsc_annot_enrichment <- function(.results = "Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output/PD_GWAS_LDSC/PD_GWAS_LDSC.results",
                                           show_plot=T,
                                           save_plot=F,
                                           title = "LDSC Heritability Enrichment",
                                           subtitle = "PD GWAS"){
  res <- data.table::fread(.results)
  res$Category <- gsub("*_0$","", res$Category)
  res$Group <- gsub("^([^_]*_[^_]*)_.*$", "\\1", res$Category)
  # Get rid of absurd enrichment values
  res <- subset(res, Enrichment>-5000 & Enrichment<5000)
  res$Valence <- ifelse(res$Enrichment>=0,1,-1)
  res$p_adj <- p.adjust(res$Enrichment_p, method="fdr")


  POLYFUN.annot_enrichment_plot  <- function(res,
                                             title = "LDSC Heritability Enrichment",
                                             subtitle = "PD GWAS"){
    sig_res <- subset(res, p_adj<0.05)
    nudge_x <- ifelse(nrow(res)>150, -.5, -.03)
    gp <- ggplot(res, aes(x=Category, y=Enrichment, fill=Group)) +
      geom_col(show.legend = F) +
      geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error),
                    width=.5, position=position_dodge(.9)) +
      geom_text(data = sig_res, aes(x=Category, y=(Enrichment+Enrichment_std_error*Valence)+5*Valence),
                label="*", nudge_x=nudge_x, color="magenta", size=7) +
      coord_flip() +
      labs(title = title,
           subtitle = subtitle) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    if(nrow(res)>150){ gp <- gp + theme(axis.text.y=element_text(size=7)) }
    return(gp)
  }

  if(show_plot){
    hist(res$Enrichment, breaks = 50)
    hist(res$Enrichment_p, breaks = 50)
    plot_dir <- dirname(.results)
    # Alphabetical
    gp.all <- POLYFUN.annot_enrichment_plot(res, title, subtitle)
    if(show_plot){print(gp.all)}
    if(save_plot){ggsave(gp.all, filename = file.path(plot_dir,"ldsc_annot_enrich_all.png"),
                         dpi = 400, width=10, height=20)}
    # Top p-vals
    gp.sig <- POLYFUN.annot_enrichment_plot(subset(res, p_adj<0.05), title, subtitle)
    if(show_plot){print(gp.sig)}
    if(save_plot){ggsave(gp.sig, filename = file.path(plot_dir,"ldsc_annot_enrich_sig.png"),
                         dpi = 400, width=10,height=5)}
  }


  POLYFUN.get_annot_refs <- function(res,
                                     supp_file="./echolocatoR/tools/polyfun/SuppTables.xlsx", sheet="S1"){
    supp <- openxlsx::read.xlsx(supp_file,sheet = sheet)
    supp_merge <- data.table:::merge.data.table(data.table::data.table(supp), res,
                                  by.x = "Annotation", by.y = "Category")
    return(supp_merge)
  }
  # supp_merge <- POLYFUN.get_annot_refs(sig_res)
  # createDT(supp_merge[,c("Annotation","Reference","Prop._SNPs","Prop._h2")])
 return(res)
}




#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.gather_annot_proportions <- function(base_url="/sc/arion/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB"){
  prop.file <- file.path(base_url,"annot_proportions.csv")
  if(!file.exists(prop.file)){
    all_annot <- list.files(base_url, pattern = "*.annot.parquet$", full.names = T)
  annot_PROP <- parallel::mclapply(all_annot, function(x){
    print(x)
    annot <- POLYFUN.read_parquet(parquet_path = x)
    annot_sums <- colSums(annot[,6:ncol(annot)])
    annot_prop <- annot_sums/nrow(annot)
    return(annot_prop)
  }, mc.cores = 4)

  annot_PROP_DF <- do.call("cbind",annot_PROP) %>% `colnames<-`(paste0("chr",1:length(annot_PROP)))
  write.csv(annot_PROP_DF, file = prop.file)
  } else {
    annot_PROP_DF <- read.csv(prop.file, row.names = 1)
  }
  return(annot_PROP_DF)
}



#' Run functional enrichment tests
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.functional_enrichment <- function(finemap_dat,
                                          PP_thresh=.95,
                                          save_plot="./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/annot_enrichment.png"){
  # "...functional enrichment of fine-mapped common SNPs in the PIP range,
  ## defined as the proportion of common SNPs in the PIP range lying in the annotation
  ## divided by the proportion of genome-wide common SNPs lying in the annotation"
  base_url <- "/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB"
  chrom <- finemap_dat$CHR[1]
  annot.file <- list.files(base_url, pattern = paste0(".UKB.",chrom,".annot.parquet"), full.names = T)

  library(reticulate)
  CONDA.activate_env(conda_env = "echoR")
  # pd <- reticulate::import("pandas")
  annot <- POLYFUN.read_parquet(parquet_path = annot.file,
                                conda_env = "echoR",
                                method = "pandas")
  annot_names <- annot %>% dplyr::select(-c(SNP,CHR,BP,A1,A2)) %>% colnames()
  annot_DT <- data.table:::merge.data.table(finemap_dat,
                                            data.table::data.table(annot) %>%
                                              dplyr::rename(SNP_y = SNP, A1_y=A1, A2_y=A2),
                                            by.x = c("CHR","POS"),
                                            by.y = c("CHR","BP"),
                                            all.x = T)
  ## SNP Groups
  # Nominal sig. GWAS
  nom.sig.GWAS <- annot_DT %>% dplyr::summarise_at(.vars = vars(annot_names),
                                                  .funs = list(sum(P < 0.05 & .>0) / n())) %>% t() %>% `colnames<-`("nom.sig.GWAS")
  # Sig. GWAS
  sig.GWAS <- annot_DT %>% dplyr::summarise_at(.vars = vars(annot_names),
                                                .funs = list(sum(P < 5e-8 & .>0) / n())) %>% t() %>% `colnames<-`("sig.GWAS")
  # Lead GWAS
  lead.GWAS <- annot_DT %>% dplyr::summarise_at(.vars = vars(annot_names),
                                               .funs = list(sum(leadSNP==T & .>0) / n())) %>% t() %>% `colnames<-`("lead.GWAS")
  # Across all CS
  UCS <- annot_DT %>% dplyr::summarise_at(.vars = vars(annot_names),
                                          .funs = list(sum(Support > 0 & .>0) / n())) %>% t() %>% `colnames<-`("UCS")
  # Tool-specific CS
  FM_methods <- gsub("*\\.PP$","",grep(pattern = "*\\.PP$",colnames(finemap_dat), value = T))
  CS <- lapply(FM_methods, function(m){
    print(m)
    annot_mod <- annot_DT
    # col_ind <- grep(paste0("^",m,".PP$"),colnames(annot_mod))
    # colnames(annot_mod)[col_ind] <- "target.PP"
    annot_mod <- annot_mod %>% dplyr::select(paste0(m,".PP"),"Support",annot_names)
    colnames(annot_mod)[1] <- "PP"
    annot_mod %>% dplyr::summarise_at(.vars = vars(annot_names),
                                     .funs = list( sum(PP > PP_thresh & .>0) / n() ))
  }) %>% data.table::rbindlist() %>% t() %>% `colnames<-`(FM_methods)

  # Gather annotation proportions genome-wide
  annot_SUMS_DF <- POLYFUN.gather_annot_proportions(base_url)
  SNP_groups <- cbind(nom.sig.GWAS, sig.GWAS, lead.GWAS, UCS, CS)
  enrich <- (SNP_groups / rowMeans(annot_SUMS_DF)) %>% melt() %>%
    `colnames<-`(c("Annot","SNP_group","Enrichment"))

  # Plot
  gp <- ggplot(enrich, aes(x=Annot, y=Enrichment, fill=SNP_group)) +
    geom_col(position="dodge", show.legend = F) + coord_flip() +
    facet_grid(facets = . ~ SNP_group) + #scales = "free_x"
    theme_bw()+ theme(axis.text.y = element_text(size = 7))
  print(gp)
  if(save_plot!=F){
    ggsave(filename = save_plot, plot = gp, dpi = 400, height = 20, width = 12)
  }
  return(enrich)
}





#' Run heritability enrichment tests
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
POLYFUN.h2_enrichment <- function(h2_df,
                                  target_SNPs=NULL,
                                  fillNA=T){
  # Only consider SNPs that overlap between LDCS and GWAS results to make things fair
  # target_SNPs <- intersect(target_SNPs, h2_df$SNP)
  if(fillNA) h2_df[is.na(h2_df$SNPVAR),"SNPVAR"] <- 0
  h2.target <- subset(h2_df, SNP %in% target_SNPs)
  if(nrow(h2.target)>0){
    # Calculate enrichment
    target_h2 <- sum(h2.target$SNPVAR, na.rm = T)
    total_h2 <- sum(h2_df$SNPVAR, na.rm = T)
    n_target_SNPs <- nrow(h2.target)
    n_total_SNPs <- nrow(h2_df)

    h2.enrichment <- (target_h2/total_h2) / (n_target_SNPs/n_total_SNPs)
  } else {
    h2.enrichment <- NA
  }
  return(h2.enrichment)
}
# POLYFUN.h2_enrichment <- function(h2_df,
#                                   target_SNPs=NULL){
#   h2.target <- subset(h2_df, SNP %in% unique(target_SNPs))
#   data.table:::merge.data.table(subset(finemap_dat))
# }





#' Run heritability enrichment tests across SNP groups
#' @source
#' https://www.biorxiv.org/content/10.1101/807792v3
#' @keywords internal
#' @family polyfun
#' @examples
#' \dontrun{
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide"
#' # IMPORTANT! For this to make sense, you need to merge the full data ("merged_DT" only includes Support>0 and leadSNPs)
#' merged_dat <- merge_finemapping_results(dataset = dirname(root), LD_reference = "UKB", minimum_support = 0)
#' merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)
#'
#' RES <- POLYFUN.h2_enrichment_SNPgroups(merged_dat=merged_dat, ldsc_dir=file.path(root,"PolyFun/output"),  save_enrich=file.path(root,"PolyFun/Nalls23andMe_2019.h2_enrich.snp_groups.csv.gz"))
#' }
POLYFUN.h2_enrichment_SNPgroups <- function(merged_dat,
                                            chrom="*",
                                            ldsc_dir="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output",
                                            ldsc_suffix="*.snpvar_constrained.gz",
                                            save_enrich=F,
                                            nThread=4){
  # Gather your heritability
  ldsc.files <- list.files(ldsc_dir, pattern = ldsc_suffix, full.names = T) %>%
    grep(pattern = paste0(".",chrom,"."), value = T)
  h2_DF <- echolocatoR:::.rbind.file.list(ldsc.files)

  h2_merged <- data.table::merge.data.table(merged_dat,
                                            subset(h2_DF, select=c(SNP,SNPVAR)),
                                            all.x = T,
                                            by=c("SNP"))
  if(!"Support_noPF" %in% colnames(merged_dat)){
    merged_dat <- find_consensus_SNPs_no_PolyFun(merged_dat)
  }

  # Iterate over loci
  RES <-parallel::mclapply(unique(merged_dat$Locus), function(locus){
    print(locus)
    finemap_dat <- subset(merged_dat, Locus==locus)
    h2_df <- subset(h2_DF, SNP %in% unique(finemap_dat$SNP))

    # Random subset  of the same size as the Consenus SNPs
    size <- dplyr::n_distinct(subset(finemap_dat, Consensus_SNP)$SNP)
    random <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                    target_SNPs= sample(finemap_dat$SNP, size = size) )
    # All SNPs
    GWAS.all <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs=finemap_dat$SNP )
    # GWAS nominally sig hits
    GWAS.nom.sig <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                          target_SNPs=subset(finemap_dat, P<.05)$SNP )
    # GWAS sig hits
    GWAS.sig <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs=subset(finemap_dat, P<5e-8)$SNP)
    # Lead GWAS
    GWAS.lead <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                       target_SNPs=subset(finemap_dat, leadSNP)$SNP)
    # Credible Set
    UCS <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                 target_SNPs = subset(finemap_dat, Support>0)$SNP)
    # # PAINTOR CS
    # PAINTOR.credset <- POLYFUN.h2_enrichment(h2_df=h2_df,
    #                                          target_SNPs = subset(finemap_dat, PAINTOR.CS>0)$SNP)
    ABF.credset <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                         target_SNPs = subset(finemap_dat, ABF.CS>0)$SNP)

    FINEMAP.credset <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                             target_SNPs = subset(finemap_dat, FINEMAP.CS>0)$SNP)

    SUSIE.credset <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                           target_SNPs = subset(finemap_dat, SUSIE.CS>0)$SNP)

    POLYFUN_SUSIE.credset <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                                   target_SNPs = subset(finemap_dat, POLYFUN_SUSIE.CS>0)$SNP)
    # Support levels
    ## Support==0
    support0 <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs = subset(finemap_dat, Support==0)$SNP)
    support1 <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs = subset(finemap_dat, Support==1)$SNP)
    support2 <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs = subset(finemap_dat, Support==2)$SNP)
    support3 <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs = subset(finemap_dat, Support==3)$SNP)
    support4 <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs = subset(finemap_dat, Support==4)$SNP)

    # Consenus SNPs
    Finemap.consensus <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                               target_SNPs = subset(finemap_dat, Consensus_SNP)$SNP)

    # Consensus SNPs (no PolyFun)
    Finemap.consensus_noPF <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                                    target_SNPs = subset(finemap_dat, Consensus_SNP_noPF)$SNP)
    # UCS SNPs (no PolyFun)
    UCS_noPF <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                      target_SNPs = subset(finemap_dat, Support_noPF>0)$SNP)

    res <- data.frame(SNP_group=c("Random",
                                  "All",
                                  "GWAS nom. sig.",
                                  "GWAS sig.",
                                  "GWAS lead",
                                  "ABF CS",
                                  "SUSIE CS",
                                  "POLYFUN-SUSIE CS",
                                  "FINEMAP CS",
                                  "UCS (-PolyFun)",
                                  "UCS",
                                  "Support==0",
                                  "Support==1",
                                  "Support==2",
                                  "Support==3",
                                  "Support==4",
                                  "Consensus (-PolyFun)",
                                  "Consensus"),
                      h2.enrichment=c(random,
                                      GWAS.all,
                                      GWAS.nom.sig,
                                      GWAS.sig,
                                      GWAS.lead,
                                      ABF.credset,
                                      SUSIE.credset,
                                      POLYFUN_SUSIE.credset,
                                      FINEMAP.credset,
                                      UCS_noPF,
                                      UCS,
                                      support0,
                                      support1,
                                      support2,
                                      support3,
                                      support4,
                                      Finemap.consensus_noPF,
                                      Finemap.consensus))
    res <- cbind(Locus=locus, res)
    return(res)
  }, mc.cores = nThread) %>% data.table::rbindlist(fill = T)

  if(save_enrich!=F){
    printer("POLFUN:: Saving enrichment results ==>",save_enrich)
    dir.create(dirname(save_enrich), showWarnings = F, recursive = T)
    data.table::fwrite(RES, save_enrich, nThread=nThread)
  }
  return(RES)
}




#' Plot heritability (h2) enrichment
#'
#' @family polyfun
#' @examples
#' \dontrun{
#' root <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping"
#' merged_dat <- merge_finemapping_results(dataset = file.path("Data/GWAS/Nalls23andMe_2019"), LD_reference = "UKB", minimum_support = 0)
#' RES <- POLYFUN.h2_enrichment_SNPgroups(merged_dat=merged_dat, out.path="/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output")
#'
#' plot.h2 <- POLYFUN.h2_enrichment_SNPgroups_plot(RES = RES, show_plot = T)
#' }
POLYFUN.h2_enrichment_SNPgroups_plot <- function(RES,
                                                 snp_groups=c("GWAS lead","UCS","Consensus (-PolyFun)","Consensus"),
                                                 comparisons_filter=function(x){if("Consensus" %in% x) return(x)},
                                                 method="wilcox.test",
                                                 remove_outliers=T,
                                                 title= "S-LDSC heritability enrichment",
                                                 xlabel="SNP group",
                                                 ylabel=bquote(~'h'^2~'enrichment'),
                                                 show_xtext=T,
                                                 show_padj=T,
                                                 show_signif=T,
                                                 vjust_signif=0.5,
                                                 show_plot=T,
                                                 save_path=F){
  colorDict <- snp_group_colorDict()
  if(is.null(snp_groups)) snp_groups <- names(colorDict)
  if(remove_outliers){
    # https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
    outliers <- boxplot(RES$h2.enrichment, plot=F)$out
    RES <- RES[-which(RES$h2.enrichment %in% outliers),]
  }
  plot_dat <- subset(RES, SNP_group %in% snp_groups) %>%
    dplyr::mutate(SNP_group=factor(SNP_group, levels=names(colorDict), ordered = T))
  snp.groups <- unique(plot_dat$SNP_group)
  comparisons <- utils::combn(x = as.character(snp.groups),
                              m=2,
                              FUN = comparisons_filter,
                              simplify = F) %>% purrr::compact()
  pb <- ggplot(data = plot_dat, aes(x=SNP_group, y=h2.enrichment, fill=SNP_group)) +
    geom_jitter(alpha=.1,width = .25, show.legend = F, shape=16, height=0) +
    geom_violin(alpha=.6, show.legend = F) +
    geom_boxplot(alpha=.6, color="black", show.legend = F) +
    geom_hline(yintercept = log10(1), linetype=2, alpha=.5) +
    # scale_fill_manual(values =  c("red","green3","goldenrod3")) +
    labs(y=ylabel, x=xlabel,
         title=title) +
    theme(legend.position = "none")  +
    scale_fill_manual(values = colorDict) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  if(show_padj){
    pb <- pb + ggpubr::stat_compare_means(method = method,
                                 comparisons = comparisons,
                                 label = "p.adj", vjust=2, size=3)
  }
  if(show_signif){
    pb <- pb + ggpubr::stat_compare_means(method = method,
                               comparisons = comparisons,
                               label = "p.signif", size=3, vjust =  vjust_signif)
  }

  if(!show_xtext){
    pb <- pb + theme(axis.text.x = element_blank(),
                      axis.title.x = element_blank())
  }

  if(show_plot) print(pb)
  if(save_path!=F){
    # save_path <- '~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/snp_group.h2_enrichment.png'
    ggsave(save_path,
           plot = pb, dpi = 300, height=9, width=11)
  }
  return(pb)
}






