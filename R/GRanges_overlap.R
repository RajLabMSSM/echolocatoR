#' Find overlap between genomic coordinates/ranges
#'
#' @family GRanges
#' @keywords internal
GRanges_overlap <- function(dat1,
                            dat2,
                            chrom_col.1="chrom",
                            start_col.1="start",
                            end_col.1="end",
                            chrom_col.2="chrom",
                            start_col.2="start",
                            end_col.2="end",
                            return_merged=T,
                            chr_format="NCBI",
                            verbose=F){
  # dat1
  if(class(dat1)[1]=="GRanges"){
    printer("+ dat1 already in GRanges format", v=verbose)
    gr.dat1 <- dat1
  } else {
    gr.dat1 <- GenomicRanges::makeGRangesFromDataFrame(dat1,
                                                       seqnames.field = chrom_col.1,
                                                       start.field = start_col.1,
                                                       end.field = end_col.1,
                                                       ignore.strand = T,
                                                       keep.extra.columns = T)
  }
  # dat2
  if(class(dat2)[1]=="GRanges"){
    printer("+ dat2 already in GRanges format", v=verbose)
    gr.dat2 <- dat2
  } else{
    printer("+ Converting dat2 to GRanges", v=verbose)
    gr.dat2 <- GenomicRanges::makeGRangesFromDataFrame(dat2,
                                                       seqnames.field = chrom_col.2,
                                                       start.field = start_col.2,
                                                       end.field = end_col.2,
                                                       ignore.strand = T,
                                                       keep.extra.columns = T)
  }
  # Standardize seqnames format
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.dat1) <- chr_format)
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.dat2) <- chr_format)
  hits <- GenomicRanges::findOverlaps(query = gr.dat1,
                                      subject = gr.dat2)
  gr.hits <- gr.dat2[ S4Vectors::subjectHits(hits), ]
  if(return_merged){
    printer("+ Merging both GRanges.", v=verbose)
    GenomicRanges::mcols(gr.hits) <- cbind(GenomicRanges::mcols(gr.hits),
                                           GenomicRanges::mcols(gr.dat1[S4Vectors::queryHits(hits),]) )
  }

  # gr.hits <- cbind(mcols(gr.regions[ S4Vectors::subjectHits(hits), ] ),
  #                         mcols(gr.consensus[S4Vectors::queryHits(hits),]) )
  message("",nrow(GenomicRanges::mcols(gr.hits))," query SNP(s) detected with reference overlap." )
  # print(data.frame(mcols(gr.hits[,c("Name","SNP")])) )
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.hits) <- chr_format)
  return(gr.hits)
}


