snp_group_filters <- function(invert=F,
                              random_sample_size=20){
  if(is.null(random_sample_size)) random_sample_size <- 100
  snp_filters <-
    c("Random" = paste0("SNP %in% sample(sampling_df$SNP, size=",random_sample_size,")"),
      "All" = "!is.na(SNP)",
      "GWAS nom. sig."="P<0.05",
      "GWAS sig."="P<5e-8",
      "GWAS lead" = "leadSNP==T",

      "ABF CS"="ABF.CS>0",
      "SUSIE CS"="SUSIE.CS>0",
      "POLYFUN-SUSIE CS"="POLYFUN_SUSIE.CS>0",
      "FINEMAP CS"="FINEMAP.CS>0",
      "UCS (-PolyFun)"="Support_noPF>0",
      "UCS"="Support>0",

      "Support==0"="Support==0",
      "Support==1"="Support==1",
      "Support==2"="Support==2",
      "Support==3"="Support==3",
      "Support==4"="Support==4",

      "Consensus (-PolyFun)"="Consensus_SNP_noPF==T",
      "Consensus"="Consensus_SNP==T"
    )
  if(invert)  snp_filters <- setNames(names(snp_filters), unname(snp_filters))
  return(snp_filters)
}
