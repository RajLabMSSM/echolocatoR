

determine_chrom_type <- function(chrom_type=NULL,
                                 file_path=NULL,
                                 chrom_col="CHR",
                                 nThread=4,
                                 verbose=T){
  if(!is.null(chrom_type)){
    has_chr <- grepl("chr",chrom_type)
  } else {
    # Slow, but it works
    printer("Determining chrom type from file header")
    header <- get_header(large_file = file_path, colnames_only = FALSE)
    has_chr <- grepl("ch",tolower(header[[chrom_col]][1]))
  }
  printer("Chromosome format =",if(has_chr) "chr1" else "1", v=verbose)
  return(has_chr)
}

