#' Check genome build
#'
#' If \code{fullSS_genome_build==NULL} and \code{munged=TRUE},
#' infers genome build (hg19 vs. hg38)
#' from summary statistics using \link[MungeSumstats]{get_genome_builds}.
#' This can only be done with summary statistics that have already been
#' munged by \link[MungeSumstats]{format_sumstats}.
#' When \code{fullSS_genome_build} is a synonym of hg19 or hg38, this function
#' simply returns a standardized version of the user-provided
#' genome build.
#' @inheritParams finemap_locus
#' @inheritParams MungeSumstats::get_genome_builds
#' @returns Character string indicating genome build.
#'
#' @export
#' @examples
#' fullSS_path <- echodata::example_fullSS()
#' build <- check_genome(fullSS_genome_build="hg19",
#'                       fullSS_path=fullSS_path)
check_genome <- function(fullSS_genome_build=NULL,
                         munged=FALSE,
                         fullSS_path=NULL,
                         sampled_snps=10000,
                         names_from_paths=TRUE,
                         dbSNP=144,
                         nThread=1,
                         verbose=TRUE){

  gbuild <- fullSS_genome_build
  if(is.null(gbuild)){
    if(munged & (!is.null(fullSS_path))){
      messager("+ Inferring genome build",v=verbose)
      requireNamespace("MungeSumstats")
      gbuild <- MungeSumstats::get_genome_builds(
        sumstats_list = fullSS_path,
        sampled_snps = sampled_snps,
        names_from_paths = names_from_paths,
        dbSNP = dbSNP,
        nThread = nThread)
      return(gbuild)
    } else {
      messager("WARNING:: fullSS_genome_build not provided.",
               "Assuming 'GRCH37'.",
               v=verbose)
      gbuild <- "GRCH37"
      return(gbuild)
    }
  } else {
    opts <- list(hg19=tolower(c("hg19","hg37","GRCh37","grch37")),
                 hg38=tolower(c("hg38","GRCh38","grch38")))
    gbuild <- if(tolower(gbuild) %in% opts$hg19) "GRCH37" else gbuild
    gbuild <- if(tolower(gbuild) %in% opts$hg38) "GRCH38" else gbuild
    return(gbuild)
  }
}
