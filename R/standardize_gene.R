#' Standardize genes
#'
#' Convert ensembl IDs to gene symbols.
#'
#' @keywords internal
#' @importFrom tidyr separate
standardize_gene <- function(dat,
                             gene_col="gene",
                             verbose=TRUE){
  requireNamespace("orthogene")

  messager("Standardizing gene name",v=verbose)
  dat <- tidyr::separate(dat,
                         col=gene_col,
                         sep = ":",
                         into = c("chr","start","end","ensembl_id"),
                         fill="left",
                         remove  = FALSE)
  genes <- gsub("\\..*","",dat$ensembl_id)
  dat$Gene <- orthogene::map_genes(genes = genes)$name
  return(dat)
}

