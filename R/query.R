#' Query handler
#'
#' Handles which query method to use
#'
#' @family query functions
#' @param subset_path Path of the resulting locus subset file.
#' @inheritParams finemap_locus
#' @inheritParams finemap_loci
#'
#' @keywords internal
#' @importFrom echotabix convert_and_query
query_handler <- function(fullSS_path,
                          locus_dir=NULL,
                          top_SNPs=NULL,
                          subset_path,
                          min_POS=NA,
                          max_POS=NA,
                          bp_distance=500000,
                          locus_col="Gene",
                          chrom_col="CHR",
                          chrom_type=NULL,
                          position_col="POS",
                          file_sep="\t",
                          query_by="coordinates",
                          conda_env="echoR",
                          force_new_subset = FALSE,
                          verbose=TRUE){
  messager("+ Query Method:",query_by,v=verbose)
  locus <- basename(locus_dir)

  if(query_by=="tabix"){
    #### Subset top_SNPs to query with ####
    topSNP_sub <- top_SNPs[top_SNPs$Locus==locus & !is.na(top_SNPs$Locus),]
    if(detect_genes(loci = locus, verbose  = FALSE)){
      topSNP_sub <- subset(topSNP_sub, Gene==names(locus))
    }
    study_dir <- get_study_dir(locus_dir)

    query_res <- echotabix::convert_and_query(
       ## Target args
       target_path = fullSS_path,
       target_chrom_col = chrom_col,
       target_start_col = position_col,
       ## Query args
       query_save = TRUE,
       query_save_path = subset_path,
       study_dir = study_dir,
       query_granges = topSNP_sub,
       verbose = verbose,
       query_force_new = force_new_subset
    )
  }
  #### Munge the full sumary stats file ####
  # Read in a standardize the entire full summary stats file at once.
  if(query_by=="fullSS"){
    file.copy(fullSS_path, subset_path)
  }
  return(query_res)
}


