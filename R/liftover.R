
#' Lift genome across builds
#'
#' @param build_conversion "hg19.to.hg38" (\emph{default}) or "hg38.to.hg19.
#' @family utils
#' @examples
#' data("BST1")
#' gr.lifted <- LIFTOVER(dat=BST1, build.conversion="hg19.to.hg38")
LIFTOVER <- function(dat,
                     build.conversion="hg19.to.hg38",
                     chrom_col="CHR",
                     start_col="POS",
                     end_col="POS",
                     chr_format="NCBI",
                     return_as_granges=T,
                     verbose=T){
  printer("XGR:: Lifting genome build:", build.conversion, v = verbose)
  # Save original coordinates and SNP IDs
  dat <- dat %>% dplyr::mutate(chrom=paste0("chr",gsub("chr","",eval(parse(text=chrom_col)))) )
  dat[,paste0("POS.",strsplit(build.conversion, "\\.")[[1]][1])] <- dat[[start_col]]
  # chain <- rtracklayer::import.chain(con = chain_paths$hg19_to_hg38)
  gr.dat <- GenomicRanges::makeGRangesFromDataFrame(df = dat,
                                                    keep.extra.columns = T,
                                                    # "chrom" col is created above
                                                    seqnames.field = "chrom",
                                                    start.field = start_col,
                                                    end.field = end_col)
  if("XGR" %in% rownames( installed.packages())){
    liftover_func <- XGR::xLiftOver()
  } else {
    liftover_func <- echolocatoR::xLiftOver()
  }
  gr.lifted <- liftover_func(data.file = gr.dat,
                             format.file = "GRanges",
                             build.conversion = build.conversion,
                             verbose = F ,
                             merged = F)  # merge must =F in order to work
  # Standardize seqnames format
  GenomeInfoDb::seqlevelsStyle(gr.lifted) <- chr_format
  # Convert  to df
  if(return_as_granges==F){
    gr.lifted <- data.frame(gr.lifted) %>%
      dplyr::mutate(POS=start) %>%
      dplyr::select(-c("seqnames","start","end","width","strand"))
  }
  return(gr.lifted)
}


