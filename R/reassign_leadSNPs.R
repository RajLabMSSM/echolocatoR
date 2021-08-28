reassign_leadSNPs <- function(merged_dat,
                              snp_col="SNP",
                              pval_col="P",
                              effect_col="Effect",
                              grouping_vars=c("Dataset","Locus"),
                              nThread=1,
                              verbose=T){
  # dplyr::mutate(merged_dat, id=paste(eval(parse(text=)), sep="."))
  merged_dat$id <- paste(merged_dat$Dataset,merged_dat$Locus,sep=".")
  merged_tmp <- dplyr::rename(merged_dat,
                              SNP=snp_col,
                              P=pval_col,
                              Effect=effect_col)
  merged_dat$leadSNP <- F
  MERGED_DAT <- parallel::mclapply(unique(merged_dat$id), function(ID){
    printer(ID, v=verbose)
    lead.snps <- (merged_tmp[merged_tmp$id==ID,] %>%
                    # dplyr::group_by(Dataset, Locus) %>%
                    dplyr::arrange(P, desc(abs(Effect))) %>%
                    dplyr::slice(1))$SNP %>% unique()
    merged_dat[merged_dat$id==ID,"leadSNP"] <- merged_dat[merged_dat$id==ID,][[snp_col]] %in% lead.snps
    return(merged_dat)
  }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)

  if(!sum(MERGED_DAT$leadSNP)==dplyr::n_distinct(MERGED_DAT$id)) warning("leadSNP count doesn't match up with unique id count.", v=verbose)
  return(MERGED_DAT)
}



