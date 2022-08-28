#' Extract a subset of the summary stats
#'
#' Use \emph{tabix} to extract a locus subset
#'  from the full summary statistics file.
#'
#' @inheritParams finemap_locus
#' @inheritParams echodata::standardize
#'
#' @family query functions
#' @keywords internal
#' @importFrom echodata check_if_empty standardize
extract_snp_subset <- function(subset_path,
                               locus = NULL,
                               colmap=echodata::construct_colmap(),
                               fullSS_path,
                               topSNPs,
                               LD_reference,
                               force_new_subset=FALSE,
                               bp_distance=500000,
                               superpopulation="EUR",
                               compute_n="ldsc",

                               query_by="tabix",
                               nThread=1,
                               conda_env="echoR_mini",
                               verbose=TRUE){

  # echoverseTemplate:::source_all();
  # echoverseTemplate:::args2vars(extract_snp_subset);

  #### Set up paths ####
  locus_dir <- get_locus_dir(subset_path = subset_path)
  if(is.null(locus)){
    locus <- basename(locus_dir)
  }
  multi_path <- echofinemap::create_method_path(
    locus_dir = locus_dir,
    LD_reference = LD_reference,
    finemap_method = "Multi-finemap",
    compress = TRUE)

  #### Priority 1: Check if multi-finemap results exists ####
  if (file.exists(multi_path) & isFALSE(force_new_subset)){
    messager("+ Importing pre-existing file:",multi_path, v=verbose)
    echodata::check_if_empty(multi_path)
    query <- data.table::fread(multi_path)
    return(query)

    #### Priority 2: Check is subset exists ####
  } else if(file.exists(subset_path) & isFALSE(force_new_subset)){
      messager("+ Importing pre-existing file:",subset_path, v=verbose)
      echodata::check_if_empty(subset_path)
      query <- data.table::fread(subset_path)
      return(query)

  } else {
    #### Priority 3: Convert and query ####
    query <- query_handler(locus_dir=locus_dir,
                           locus = locus,
                           fullSS_path=fullSS_path,
                           subset_path=subset_path,
                           topSNPs=topSNPs,
                           colmap=colmap,
                           bp_distance=bp_distance,
                           query_by=query_by,
                           force_new_subset=force_new_subset,
                           nThread=nThread,
                           conda_env=conda_env,
                           verbose=verbose)
    #### Standardize ####
    query <- echodata::standardize(query=query,
                                   subset_path=subset_path,
                                   colmap=colmap,
                                   locus=locus,
                                   compute_n=compute_n,
                                   nThread=nThread,
                                   verbose=verbose)
  }
  #### Return ####
  return(query)
}
