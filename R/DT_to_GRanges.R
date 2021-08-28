#' Convert data.table to GRanges object
#'
#' @family XGR
#' @keywords internal
DT_to_GRanges <- function(subset_DT,
                          style="NCBI",
                          chrom_col="CHR",
                          position_col="POS"){
  subset_DT[["SEQnames"]] <- subset_DT[[chrom_col]]
  gr.snp <- biovizBase::transformDfToGr(subset_DT,
                                        seqnames = "SEQnames",
                                        start = "POS",
                                        end = "POS")
 suppressMessages(suppressWarnings( GenomeInfoDb::seqlevelsStyle(gr.snp) <- style))
  return(gr.snp)
}
