
###### Tutorial
# https://bioconductor.org/packages/release/bioc/vignettes/cpvSNP/inst/doc/cpvSNP.pdf

FM_all <- merge_finemapping_results(minimum_support = 0, dataset = "./Data/GWAS/Nalls23andMe_2019")
FM_gene <- subset(FM_all, Gene=="LRRK2", select=c("P","SNP","POS","CHR"))
results_path <- "./Data/GWAS/Nalls23andMe_2019/LRRK2/"
LD_matrix <- readRDS(file.path(results_path, "plink/LD_matrix.RData"))

library(cpvSNP) # BiocManager::install("cpvSNP")
data(geneSetAnalysis)
names(geneSetAnalysis) 
arrayDataGR <- createArrayData(FM_gene, 
                               positionName="POS",
                               chromosomeName = "CHR")
# geneSets <- geneSetAnalysis[["geneSets"]]
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genesHg19 <- genes(txdb)
geneSets <- geneSetAnalysis[["geneSets"]]

set1 <- GSEABase::GeneSet("LRRK2",geneIdType=SymbolIdentifier(), setName="LRRK2_only")
set2 <- GSEABase::GeneSet(unique(FM_all$Gene),geneIdType=SymbolIdentifier(), setName="all_PD_loci")

gene_sets <- GSEABase::GeneSetCollection(set1, set2)

snpsGSC <- geneToSNPList(geneList = geneSets, 
                         arrayData = arrayDataGR, 
                         genes = genesHg19)


# ldMat <- geneSetAnalysis[["ldMat"]]

vRes <- vegas(snpsGSC[1], arrayDataGR, LD_matrix)
vRes
summary(unlist(simulatedStats(vRes)))
pValue(vRes)
degreesOfFreedom(vRes)
statistic(vRes)
