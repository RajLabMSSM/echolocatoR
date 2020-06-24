## ---- include = FALSE---------------------------------------------------------
root.dir <- "~/Desktop"
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  root.dir = root.dir,
  fig.height = 12,
  fig.width = 10
)  
knitr::opts_knit$set(root.dir = root.dir)
# knitr::opts_chunk$get("root.dir")

## Build vignettes
## This will generate a /doc folder with the lmotted html files.
## Go to GitHub and enable GitHub pages, selecting "branch/doc"
# devtools::build_vignettes(quiet = F, clean=F)
# https://resources.github.com/whitepapers/github-and-rstudio/

# devtools::build()

## ----setup, root.dir="~/Desktop"----------------------------------------------
library(echolocatoR) 

## ----Prepare `top_SNPs` data.frame--------------------------------------------
data("Nalls_top_SNPs");
top_SNPs <- import_topSNPs(
  # topSS = "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx",
  topSS = Nalls_top_SNPs,
  chrom_col = "CHR", position_col = "BP", snp_col="SNP",
  pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
  group_by_locus = T,
  locus_col = "Nearest Gene",
  remove_variants = "rs34637584") 

## ----fullSS-------------------------------------------------------------------
fullSS_path <- example_fullSS(fullSS_path="~/Desktop/Nalls23andMe_2019.fullSS_subset.tsv")

## ----Run fine-mapping pipeline------------------------------------------------
Nalls23andMe_2019.results <- finemap_loci(# GENERAL ARGUMENTS 
                                          top_SNPs = top_SNPs,  
                                          results_dir = "~/Desktop/results",
                                          loci = c("BST1","MEX3C"),#top_SNPs$Gene, 
                                          dataset_name = "Nalls23andMe_2019",
                                          dataset_type = "GWAS",  
                                          force_new_subset = T,
                                          force_new_LD = F,
                                          force_new_finemap = F,
                                          remove_tmps = T,
                                          
                 # SUMMARY STATS ARGUMENTS
                 fullSS_path = fullSS_path,
                 query_by ="tabix",
                 chrom_col = "CHR", position_col = "POS", snp_col = "RSID",
                 pval_col = "p", effect_col = "beta", stderr_col = "se",
                 freq_col = "freq", MAF_col = "calculate",
                 A1_col = "A1",
                 A2_col = "A2",
                 
                 # FILTERING ARGUMENTS
                 bp_distance = 500000*2,
                 min_MAF = 0.001, 
                 trim_gene_limits = F,
                 
                 # FINE-MAPPING ARGUMENTS
                 finemap_methods = c("ABF","SUSIE","FINEMAP"),
                 n_causal = 5,
                 PP_threshold = .95,
                 
                 # LD ARGUMENTS 
                 LD_reference = "UKB",#"1KG_Phase1",
                 superpopulation = "EUR",
                 LD_download_method = "axel",
                 
                 # PLOT ARGUMENTS 
                 ## general   
                 plot.types=c("fancy"),
                 plot.window = 100000,
                 ## XGR
                 # plot.XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes"), 
                 ## Roadmap
                 plot.Roadmap = F,
                 plot.Roadmap_query = NULL,
                 # Nott et al. (2019)
                 plot.Nott_epigenome = T, 
                 plot.Nott_binwidth = 100
                 )

