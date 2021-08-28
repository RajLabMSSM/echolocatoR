#' Ensembl IDs -> Gene symbols
#'
#' Convert HGNC gene symbols to Ensembl IDs.
#'
#' @param gene_symbols List of HGNC gene symbols.
#' @return List of ensembl IDs.
#' @examples
#' gene_symbols <- c("FOXP2","BDNF","DCX","GFAP")
#' ensembl_IDs <- hgnc_to_ensembl(gene_symbols)
#' @keywords internal
hgnc_to_ensembl <- function(gene_symbols){
  gene_symbols[is.na(gene_symbols)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = gene_symbols,
                                      keytype = "SYMBOL",
                                      column = "GENEID")
  return(conversion)
}




#' Gene symbols -> Ensembl IDs
#'
#' Convert Ensembl IDs to HGNC gene symbols.
#'
#' @param ensembl_ids List of ensembl IDs.
#' @return List of HGNC gene symbols.
#' @examples
#' ensembl_IDs <- c("ENSG00000128573","ENSG00000176697","ENSG00000077279","ENSG00000131095" )
#' gene_symbols <- ensembl_to_hgnc(ensembl_IDs)
#' @keywords internal
ensembl_to_hgnc <- function(ensembl_ids){
  ensembl_ids <- gsub("\\..*","",ensembl_ids) # Remove transcript suffixes
  ensembl_ids[is.na(ensembl_ids)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = ensembl_ids,
                                      keytype = "GENEID",
                                      column = "SYMBOL")
  return(conversion)
}


