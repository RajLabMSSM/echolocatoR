




LD.compare_panels <- function(root="/sc/arion/projects/pd-omics/brian/Fine_Mapping",
                              dataset="Data/GWAS/Nalls23andMe_2019"){
  # merged_UKB <- data.table::fread(file.path(root,"Data/GWAS/Nalls23andMe_2019/_genome_wide/merged_UKB.csv.gz"))
  merged_UKB <- merge_finemapping_results(dataset = file.path(root,dataset),
                                          LD_reference = "UKB",
                                          minimum_support = 1, include_leadSNPs = T)
  # merged_1KG <- data.table::fread(file.path(root,"Data/GWAS/Nalls23andMe_2019/_genome_wide/merged_1KGphase3.csv.gz"))
  merged_1KG <- merge_finemapping_results(dataset = dataset,
                                          LD_reference = "1KGphase3",
                                          minimum_support = 1, include_leadSNPs = T)
  # Merge and remove non-overlapping loci
  common_loci <- dplyr::intersect(unique(merged_UKB$Locus),
                                  unique(merged_1KG$Locus))
  # merged_LD <- rbind(subset(merged_UKB, Locus %in% common_loci) %>%
  #                      dplyr::mutate(LD_panel="UKB"),
  #                    subset(merged_1KG, Locus %in% common_loci) %>%
  #                      dplyr::mutate(LD_panel="1KG")
  #                    )
  merged_LD <- data.table::merge.data.table(#UKB
                                            subset(merged_UKB, Locus %in% common_loci) %>%
                                               `colnames<-`(paste(colnames(.),"UKB",sep='.')) %>%
                                              dplyr::rename(Locus=Locus.UKB, SNP=SNP.UKB),
                                            # 1KG
                                            subset(merged_1KG, Locus %in% common_loci) %>%
                                               `colnames<-`(paste(colnames(.),"1KG",sep='.')) %>%
                                              dplyr::rename(Locus=Locus.1KG, SNP=SNP.1KG),
                                            by = c("Locus","SNP")
                                             )
   # overlap <- merged_LD %>%
   #   dplyr::group_by(Locus) %>%
   #   dplyr::summarise(leadSNP_overlap=sum(leadSNP.UKB>0 & leadSNP.1KG>0, na.rm = T),
   #                    leadSNP_prop=sum(leadSNP.UKB>0 & leadSNP.1KG>0, na.rm = T)/sum(leadSNP.UKB>0 | leadSNP.1KG>0, na.rm = T),
   #
   #                    UCS_overlap=sum(Support.UKB>0 & Support.1KG>0, na.rm = T),
   #                    UCS_prop=sum(Support.UKB>0 & Support.1KG>0, na.rm = T)/sum(Support.UKB>0 | Support.1KG>0, na.rm = T),
   #
   #                    Consensus_overlap=sum(Consensus_SNP.UKB>0 & Consensus_SNP.1KG>0, na.rm = T),
   #                    Consensus_prop=sum(Consensus_SNP.UKB>0 & Consensus_SNP.1KG>0, na.rm = T)/sum(Consensus_SNP.UKB>0 | Consensus_SNP.1KG>0, na.rm = T)
   #                    ) %>% data.frame()
   # colMeans(overlap[,-1], na.rm = T)


   overlap <- merged_LD %>%
     dplyr::mutate(leadSNP_overlap=leadSNP.UKB>0 & leadSNP.1KG>0,
                   UCS_overlap=Support.UKB>0 & Support.1KG>0,
                   Consensus_overlap=Consensus_SNP.UKB>0 & Consensus_SNP.1KG>0)



   overlap %>%
     dplyr::group_by(Locus) %>%
     dplyr::summarise(leadSNP_overlap=sum(leadSNP_overlap, na.rm = T),
                      leadSNP_prop=sum(leadSNP_overlap, na.rm = T)/sum(leadSNP_overlap, na.rm = T),

                      UCS_overlap=sum(UCS_overlap, na.rm = T),
                      UCS_prop=sum(UCS_overlap, na.rm = T)/sum(UCS_overlap, na.rm = T),

                      Consensus_overlap=sum(Consensus_overlap, na.rm = T),
                      Consensus_prop=sum(Consensus_overlap, na.rm = T)/sum(Consensus_overlap, na.rm = T)
     ) %>% data.frame()



}
