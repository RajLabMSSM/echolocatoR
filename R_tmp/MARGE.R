



#' Run \emph{HOMER} TFBM enrichment
#'
#' Run HOMER transcription factor binding motif enrichment analysis
#' via the R package \code{marge}
#' @family TFBM
#' @source
#' \url{https://github.com/robertamezquita/marge}
#' @examples
#' \dontrun{
#' # devtools::install_github('robertamezquita/marge', ref = 'master')
#' ## WARNING!: Must install HOMER separately first.
#' ## ml homer R/3.6.0
#' library(marge)
#' data("merged_DT")
#' finemap_dat <- subset(merged_DT, Support>0 | leadSNP)
#' MARGE.find_motifs()
#' }
MARGE.find_motifs <- function(finemap_dat,
                              output_dir="./HOMER_results",
                              genome="hg19"){
  # chromosome, start, end, region ID
  dir.create(output_dir,showWarnings = FALSE, recursive = T)
  finemap_dat <- finemap_dat %>%
    dplyr::select(chromosome="CHR",start="POS",end="POS", dplyr::everything())
  # Run
  marge::find_motifs_genome(x = finemap_dat,
                            path = output_dir,
                            genome = genome)
}



