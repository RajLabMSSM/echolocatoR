

library(dplyr)
dat <- data.table::fread("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz")

monocyte_eqtl <- dat[grep("monocyte|myeloid",dat$qtl.id),] %>%
  subset(PP.H4>.8 & !startsWith(eGene,"RP1"))

(monocyte_eqtl %>%
    group_by(Locus.GWAS) %>%
    summarise(eGenes=dplyr::n_distinct(eGene)))$eGenes %>% summary()

data.table::fwrite(monocyte_eqtl[,c("Locus.GWAS","eGene","qtl.id","PP.H4")],  "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/colocalized_monocyte_eQTL.csv",sep=",")
