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
#' @importFrom echotabix convert_and_query construct_query
query_handler <- function(fullSS_path,
                          colmap,
                          locus_dir = NULL,
                          locus = NULL,
                          topSNPs,
                          subset_path,
                          bp_distance = 500000,
                          query_by = c("tabix","fullSS"),
                          force_new_subset = FALSE,
                          conda_env = "echoR_mini",
                          nThread = 1,
                          verbose = TRUE){

  query_by <- tolower(query_by)[1]
  messager("+ Query Method:",query_by, v=verbose)
  if(is.null(locus)){
    locus <- basename(locus_dir)
  }
  #### Subset Locus ####
  topSNP_sub <- topSNPs[topSNPs$Locus==locus & !is.na(topSNPs$Locus),]
  #### Subset to only one gene ####
  if(echodata::detect_genes(loci = locus)){
    if(!"Gene" %in% colnames(topSNP_sub)){
      messager("QTL locus detected, but Gene col not present in topSNPs.",
               "Skipping Gene filtering step.",
               v=verbose)
    } else {
      messager("QTL locus detected, subsetting to data Gene:",names(locus),
               v=verbose)
      topSNP_sub <- subset(topSNP_sub, Gene==names(locus))
    }
  }
  #### Create query GRanges ####
  query_granges <- echotabix::construct_query(
    query_chrom = topSNP_sub$CHR,
    query_start_pos = topSNP_sub$POS - bp_distance,
    query_end_pos = topSNP_sub$POS + bp_distance
  )
  #### Query ####
  if(query_by=="tabix"){
    study_dir <- get_study_dir(locus_dir = locus_dir)
    #### Tabix-index (is not already so) and query subset ####
    query <- echotabix::convert_and_query(
       ## Target args
       target_path = fullSS_path,
       target_chrom_col = colmap$CHR,
       target_start_col = colmap$POS,
       ## Query args
       query_save = TRUE,
       query_save_path = subset_path,
       study_dir = study_dir,
       query_granges = query_granges,
       query_force_new = force_new_subset,
       conda_env = conda_env,
       nThread = nThread,
       verbose = verbose
    )
  ## Remove extra query column
  if("query" %in% colnames(query)) query$query <- NULL

  #### Munge the full summary stats file ####
  } else if(query_by=="fullss"){
    ## Read in a standardize the entire
    ## full summary stats file at once.
    file.copy(fullSS_path, subset_path)
  }
  #### Report ####
  messager("+ Query:",
           formatC(dim(query)[1], big.mark = ","), "SNPs",
           "x",
           formatC(dim(query)[2],big.mark = ","),
           "columns.",
           v=verbose)
  return(query)
}
