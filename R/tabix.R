
#' Covert and query
#'
#' If it is not tabix format already
#' (determined by checking for a \link{.tbi} file of the same name in the same directory),
#' the full summary statistics file is converted into tabix format for super fast querying.
#' A query is then made using the min/max genomic positions to extract a locus-specific summary stats file.
#'
#' @family query functions
#' @inheritParams finemap_pipeline
#' @return data.table of locus subset summary statistics
#' @examples
#' \dontrun{
#' data("Nalls_top_SNPs"); data("top_SNPs")
#' fullSS_path <- example_fullSS()
#' top_SNPs_BST1 <- subset(top_SNPs, Locus=='BST1')
#' bp_distance <- 1e+06
#' min_POS <- top_SNPs_BST1$POS - bp_distance
#' max_POS <- top_SNPs_BST1$POS + bp_distance
#'
#' subset_path <- file.path("BST1_Nalls23andMe_2019_subset.tsv.gz")
#' dat <- TABIX(fullSS_path=fullSS_path,
#'              subset_path=subset_path,
#'              min_POS=min_POS,
#'              max_POS=max_POS,
#'              chrom=top_SNPs_BST1$CHR)
#' }
TABIX <- function(fullSS_path,
                  study_dir=NULL,
                  subset_path=NULL,
                  is_tabix=F,
                  chrom_col="CHR",
                  chrom_type=NULL,
                  position_col="BP",
                  min_POS,
                  max_POS,
                  chrom,
                  save_subset=T,
                  nThread=1,
                  conda_env="echoR",
                  verbose=T){
  #### Check if it's already an indexed tabix file ####
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
      file.copy(fullSS_path, tabix_out, overwrite = TRUE)
      tabix_out <- fullSS_path
    } else {
      tabix_out <- TABIX.convert_file(fullSS_path = fullSS_path,
                                      chrom_col = chrom_col,
                                      position_col = position_col,
                                      verbose=verbose)
    }
  }
  #### Check chrom format ####
  cDict <- column_dictionary(file_path = tabix_out)
  has_chr <- determine_chrom_type(chrom_type=chrom_type,
                                  file_path=tabix_out,
                                  chrom_col=chrom_col,
                                  verbose=verbose)
  chrom <- if(has_chr) paste0("chr",gsub("chr","",chrom)) else gsub("chr","",chrom)
  #### Query ####
  dat <- TABIX.query(fullSS_tabix=tabix_out,
                     chrom=chrom,
                     start_pos=min_POS,
                     end_pos=max_POS,
                     verbose=verbose)
  #### Save subset ####
  if(save_subset){
    printer("++ Saving query ==>", subset_path, v=verbose)
    dir.create(dirname(subset_path), showWarnings = FALSE, recursive = FALSE)
    data.table::fwrite(dat, file = subset_path, nThread = nThread, sep="\t")
  }
  return(dat)
}



construct_tabix_path <- function(fullSS_path,
                                 study_dir=NULL){
  study_dir <- if(is.null(study_dir)) dirname(fullSS_path) else study_dir
  fullSS.gz <- gsub(".gz|.bgz",".bgz",fullSS_path)
  tabix_out <- file.path(study_dir, basename(fullSS.gz))
  return(tabix_out)
}




infer_if_tabix <- function(file_path){
  # must meet all of these conditions in order to use a pre-existing tabix files
  file.exists(file_path) &
    (endsWith(file_path,".gz")|endsWith(file_path,".bgz")) &
    file.exists(paste0(file_path,".tbi")) &
    file.size(file_path)>0
}




#' Convert summary stats file to tabix format
#'
#' @family query functions
#' @inheritParams finemap_pipeline
#' @examples
#' \dontrun{
#' fullSS_path <- example_fullSS()
#' fullSS_tabix <- TABIX.convert_file(fullSS_path=fullSS_path, chrom_col="CHR", position_col="POS")
#' }
TABIX.convert_file <- function(fullSS_path,
                               chrom_col="CHR",
                               position_col="POS",
                               verbose=T){
  printer("TABIX:: Converting full summary stats file to tabix format for fast querying...",v=verbose)
  cDict <-  column_dictionary(file_path = fullSS_path)
  # Make sure input file isn't empty
  if(file.size(fullSS_path)==0){
    printer("TABIX:: Removing empty file =",fullSS_path);
    file.remove(fullSS_path)
  }
  ### File MUST be bgzipped first
  printer("TABIX:: Ensuring file is bgzipped.", v = verbose)
  bgz_file <- Rsamtools::bgzip(file = fullSS_path, overwrite = TRUE)

  ### Tabix-index file
  printer("TABIX:: Tabix-indexing file.")
  seqminer::tabix.createIndex(bgzipFile = bgz_file,
                              sequenceColumn = cDict[[chrom_col]],
                              startColumn = cDict[[position_col]],
                              endColumn = cDict[[position_col]],
                              ## Just use the first columns name
                              metaChar = names(cDict)[1])
  return(bgz_file)
}





#' Query a tabix file
#'
#' Query by genomic coordinates.
#'
#' @family query functions
#' @examples
#' \dontrun{
#' fullSS_path <- example_fullSS()
#' fullSS_tabix <- TABIX.convert_file(fullSS_path=fullSS_path, chrom_col="CHR", position_col="BP")
#' TABIX.query(fullSS_tabix=fullSS_tabix, chrom=4, start_pos=13737637, end_pos=13837637)
#' }
TABIX.query <- function(fullSS_tabix,
                        chrom,
                        start_pos,
                        end_pos,
                        verbose=TRUE){
  coords <- paste0(chrom,":",start_pos,"-",end_pos)
  printer("TABIX:: Extracting subset of sum stats", v=verbose)
  dat <- seqminer::tabix.read.table(tabixFile = fullSS_tabix,
                                    tabixRange = coords)
  printer("+ TABIX:: Returning",paste(dim(dat),collapse=" x "),"data.table", v=verbose)
  return(dat)
}








