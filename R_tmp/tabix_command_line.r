


construct_tabix_path <- function(fullSS_path,
                                 study_dir=NULL){
  study_dir <- if(is.null(study_dir)) dirname(fullSS_path) else study_dir
  fullSS.gz <- if(endsWith(fullSS_path,".gz")) fullSS_path else  paste0(fullSS_path,".gz")
  tabix_out <- file.path(study_dir, basename(fullSS.gz))
  return(tabix_out)
}



#' Convert summary stats file to tabix format
#'
#' @family query functions
#' @importFrom echoconda find_packages
#' @inheritParams finemap_locus
#' @examples
#' \dontrun{
#' data("genome_wide_dir");
#' fullSS_path <-  "~/Desktop/nallsEtAl2019_allSamples_allVariants.mod.txt"
#' fullSS.gz <- TABIX.convert_file(fullSS_path=fullSS_path, results_path=results_path_genome_wide, chrom_col="CHR", position_col="POS")
#' }
TABIX.convert_file <- function(fullSS_path,
                               study_dir=NULL,
                               chrom_col="CHR",
                               position_col="POS",
                               conda_env="echoR",
                               verbose=T){
  printer("TABIX:: Converting full summary stats file to tabix format for fast querying...",v=verbose)
  z_grep <- if(endsWith(fullSS_path,".gz")) "zgrep" else "grep"
  cDict <-  echodata::column_dictionary(file_path = fullSS_path) # header.path
  tabix_out <- construct_tabix_path(fullSS_path = fullSS_path,
                                       study_dir = study_dir)
  # Make sure input file isn't empty
  if(file.size(fullSS_path)==0){
    printer("TABIX:: Removing empty file =",fullSS_path);
    file.remove(fullSS_path)
  }
  cmd <- paste("(",
               # Extract the header col and sort everything else
               paste0(z_grep," '",chrom_col,"' ",fullSS_path,"; ",z_grep," -v ^'",chrom_col,"' ",fullSS_path),
               paste0("| sort -k",cDict[[chrom_col]],",",cDict[[chrom_col]]),
               paste0("-k",cDict[[position_col]],",",cDict[[position_col]],"n"),
               ")",
               # Compress with bgzip
               "|",echoconda::find_packages(packages="bgzip", conda_env=conda_env),"-f",
               ">",
               tabix_out)
  printer(cmd, v=verbose)
  system(cmd)

  # Index
  TABIX.index <- function(tabix_out,
                          chrom_i=1,
                          pos_i=2,
                          skip_lines=1,
                          conda_env="echoR",
                          verbose=T){
    tabix <- echoconda::find_packages(package="tabix",
                                conda_env=conda_env,
                                verbose = verbose)
    printer("TABIX:: Indexing",v=verbose)
    cmd2 <- paste(tabix,
                  "-f", # Force overwrite of .tbi index file
                  "-S",skip_lines,#--skip-lines
                  "-s",chrom_i,
                  "-b",pos_i,
                  "-e",pos_i,
                  tabix_out)
    printer(cmd2, v=verbose)
    system(cmd2)
  }
  TABIX.index(tabix_out=tabix_out,
              chrom_i=cDict[[chrom_col]],
              pos_i=cDict[[position_col]],
              skip_lines=1,
              conda_env=conda_env,
              verbose=verbose)
  return(tabix_out)
}




#' Query a tabix file
#'
#' Query by genomic coordinates.
#'
#' @family query functions
TABIX.query <- function(fullSS.gz,
                        chrom,
                        start_pos,
                        end_pos,
                        conda_env="echoR",
                        verbose=T){
  tabix <- echoconda::find_packages(package="tabix",
                              conda_env=conda_env,
                              verbose = verbose)
  coords <- paste0(chrom,":",start_pos,"-",end_pos)
  # cmd4 <- paste("tabix -h",fullSS.gz,coords,">",subset_path)
  printer("TABIX:: Extracting subset of sum stats", v=verbose)
  tabix_cmd <- paste(tabix,"-h",fullSS.gz,coords)
  printer("+ TABIX::",tabix_cmd, v=verbose)
  dat <- data.table::fread(cmd=tabix_cmd, nThread = 1)
  printer("+ TABIX:: Returning",paste(dim(dat),collapse=" x "),"data.table", v=verbose)
  return(dat)
}





#' Covert and query
#'
#' If it is not tabix format already
#' (determined by checking for a \link{.tbi} file of the same name in the same directory),
#' the full summary statistics file is converted into tabix format for super fast querying.
#' A query is then made using the min/max genomic positions to extract a locus-specific summary stats file.
#'
#' @family query functions
#' @inheritParams finemap_locus
#' @return data.table of locus subset summary statistics
#' @examples
#' \dontrun{
#' locus_dir <- echodata::locus_dir; data("Nalls_top_SNPs")
#' fullSS_path <- "./results/GWAS/Nalls23andMe_2019/_genome_wide/nallsEtAl2019_allSamples_allVariants.mod.txt.gz"
#' top_SNPs <- import_topSNPs(topSS = Nalls_top_SNPs, chrom_col = "CHR", position_col = "BP", snp_col="SNP", pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene", group_by_locus = T,locus_col = "Nearest Gene")
#' top_SNPs_BST1 <- subset(top_SNPs, Locus=='BST1')
#' bp_distance <- 1e+06
#' min_POS <- top_SNPs_BST1$POS - bp_distance
#' max_POS <- top_SNPs_BST1$POS + bp_distance
#' subset_file <- file.path(results_path,"BST1_Nalls23andMe_2019_subset.tsv.gz")
#' dat <- TABIX(fullSS_path=fullSS_path, subset_path="auto", min_POS=min_POS, max_POS=max_POS, chrom=top_SNPs_BST1$CHR)
#' }
TABIX <- function(fullSS_path,
                  study_dir=NULL,
                  subset_path=NULL,
                  is_tabix=FALSE,
                  chrom_col="CHR",
                  chrom_type=NULL,
                  position_col="POS",
                  min_POS,
                  max_POS,
                  chrom,
                  save_subset=T,
                  nThread=1,
                  conda_env="echoR",
                  verbose=T){
  # Check if it's already an indexed tabix file
  # fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), fullSS_path, paste0(fullSS_path,".gz"))
  tabix_out <- construct_tabix_path(fullSS_path = fullSS_path,
                                    study_dir = study_dir)
  if(infer_if_tabix(tabix_out)){
    # Checks if the file (in the study dir) already exists,
    # and whether it is a tabix-indexed file.
    printer("TABIX:: Using existing tabix file:",tabix_out,v=verbose)
    # Jump ahead and query tabix_out file
  } else {
    if(infer_if_tabix(fullSS_path)){
      printer("TABIX:: Copying existing tabix file ==>",fullSS_path,v=verbose)
      file.copy(fullSS_path, tabix_out, overwrite = T)
      tabix_out <- fullSS_path
    } else {
      tabix_out <- TABIX.convert_file(fullSS_path = fullSS_path,
                                      study_dir = study_dir,
                                      chrom_col = chrom_col,
                                      position_col = position_col,
                                      conda_env=conda_env,
                                      verbose=verbose)
    }
  }
  # Query
  cDict <- echodata::column_dictionary(file_path = tabix_out)
  # Determine chromosome format
  has_chr <- determine_chrom_type(chrom_type=chrom_type,
                                  file_path=tabix_out,
                                  chrom_col=chrom_col,
                                  nThread=nThread,
                                  verbose=verbose)
  chrom <- if(has_chr) paste0("chr",gsub("chr","",chrom)) else gsub("chr","",chrom)

  dat <- TABIX.query(fullSS.gz=tabix_out,
                     chrom=chrom,
                     start_pos=min_POS,
                     end_pos=max_POS,
                     conda_env=conda_env,
                     verbose=verbose)
  colnames(dat) <- names(cDict)
  if(save_subset){
    printer("++ Saving query ==>", subset_path, v=verbose)
    data.table::fwrite(dat, file = subset_path, nThread = nThread, sep="\t")
  }
  return(dat)
}



#
# TABIX.seqminer <- function(){
#   sure <- data.table::fread("~/Desktop/MPRA/SURE/SuRE_SNP_table_LP190708.txt.gz")
#   sure <- sure %>% dplyr::arrange(chr, SNPabspos)
#   ## This allows tabix to recognize the header
#   colnames(sure)[1] <- paste0("#",colnames(sure)[1])
#   data.table::fwrite(sure, "~/Desktop/MPRA/SURE/SuRE_SNP_table_LP190708.tsv", sep = "\t")
#
#   tab <- seqminer::tabix.createIndex("~/Desktop/MPRA/SURE/SuRE_SNP_table_LP190708.tsv.gz",
#                                      sequenceColumn = 1, startColumn = 3, endColumn = 3, metaChar = "#")
#   tab <- seqminer::tabix.read.table(tabixFile = "~/Desktop/MPRA/SURE/SuRE_SNP_table_LP190708.tsv.gz", tabixRange = "chr1:13000-20000")
# }


