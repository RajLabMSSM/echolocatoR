#' Run heritability enrichment tests
#' @source
#' https://www.nature.com/articles/s41588-020-00735-5
#' @keywords internal
#' @family polyfun
POLYFUN_h2_enrichment <- function(h2_df,
                                  target_SNPs=NULL,
                                  fillNA=0){
    # Only consider SNPs that overlap between LDCS and GWAS results to make things fair
    # target_SNPs <- intersect(target_SNPs, h2_df$SNP)
    if(!is.null(fillNA)) {
        h2_df[is.na(h2_df$SNPVAR),"SNPVAR"] <- fillNA
    }
    h2.target <- subset(h2_df, SNP %in% target_SNPs)
    if(nrow(h2.target)>0){
        # Calculate enrichment
        target_h2 <- sum(h2.target$SNPVAR, na.rm = TRUE)
        total_h2 <- sum(h2_df$SNPVAR, na.rm = TRUE)
        n_target_SNPs <- nrow(h2.target)
        n_total_SNPs <- nrow(h2_df)
        
        h2.enrichment <- (target_h2/total_h2) / (n_target_SNPs/n_total_SNPs)
    } else {
        h2.enrichment <- NA
    }
    return(h2.enrichment)
}
