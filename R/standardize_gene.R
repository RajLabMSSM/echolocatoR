standardize_gene <- function(merged_COLOC,
                             gene_col="gene",
                             verbose=T){
  printer("Standardizing gene name",v=verbose)
  merged_COLOC <- tidyr::separate(merged_COLOC,
                                  col=gene_col, sep = ":",
                                  into = c("chr","start","end","ensembl_id"),
                                  fill="left", remove = F)
  merged_COLOC$Gene <- ensembl_to_hgnc(ensembl_ids = gsub("\\..*","",merged_COLOC$ensembl_id))
  return(merged_COLOC)
}

