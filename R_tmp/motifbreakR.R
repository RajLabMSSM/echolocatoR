# motifbreakR
# https://bioconductor.org/packages/release/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html

# BiocManager::install(c("motifbreakR",
#                        "SNPlocs.Hsapiens.dbSNP142.GRCh37",
#                        "BSgenome.Hsapiens.UCSC.hg19"))
require(motifbreakR)
# Databases
library(SNPlocs.Hsapiens.dbSNP142.GRCh37) # dbSNP142 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)     # hg19 genome
# library(BSgenome)
# available.SNPs()


finemap_DT <- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt")
CS <- subset(finemap_DT, Support>0)$SNP

variants <- motifbreakR::snps.from.rsid(rsid = CS,
                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19)

motifbreakr.results <- motifbreakR::motifbreakR(snpList = variants, 
                                   pwmList = MotifDb, 
                                   threshold = 0.9, 
                                   #This can be (although not always) a very memory 
                                   # and time intensive process if the algorithm doesn’t converge rapidly.
                                   # filterp = T,
                                   verbose = T)

human_ids <- grep("HUMAN", motifbreakr.results$providerId, value = T )
human.results <- subset(motifbreakr.results, providerId %in% human_ids)
human.results <- motifbreakr.results[motifbreakr.results$providerId %in% human_ids]

motifbreakR::plotMB(results = human.results, 
                     rsid = "rs7294619", 
                     effect = "strong")
               

# rsid <- human.results[names(human.results) %in% "rs7294619"]
# calculatePvalue(rsid) 