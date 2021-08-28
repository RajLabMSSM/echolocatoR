#' Remove all SNPs outside of of a given gene.
#'
#' Get the min/max coordinates of a given gene (including known regulatory regions, introns, and exons).
#' Remove any SNPs from the data.frame that fall outside these coordinates.
#'
#' @family SNP filters
#' @keywords internal
gene_trimmer <- function(subset_DT,
                         gene,
                         min_POS=NULL,
                         max_POS=NULL){
  printer("BIOMART:: Trimming data to only include SNPs within gene coordinates.")
  gene_info <- biomart_geneInfo(gene)
  gene_info_sub <- subset(gene_info, hgnc_symbol==gene)
  # Take most limiting min position
  min_POS <- max(min_POS, gene_info_sub$start_position, na.rm = T)
  # Take most limiting max position
  max_POS <- min(max_POS, gene_info_sub$end_position, na.rm = T)
  subset_DT <- subset(subset_DT, CHR==gene_info$chromosome_name[1] & POS>=min_POS & POS<=max_POS)
  printer("BIOMART::",nrow(subset_DT),"SNPs left after trimming.")
  return(subset_DT)
}

