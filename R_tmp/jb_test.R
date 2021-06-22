library(echolocatoR)
stopifnot(packageVersion("echolocatoR") == "0.1.2")

data("Nalls_top_SNPs")
top_SNPs <- import_topSNPs(
  topSS = Nalls_top_SNPs,
  position_col = "BP",
  pval_col = "P, all studies",
  effect_col = "Beta, all studies",
  gene_col = "Nearest Gene",
  locus_col = "Nearest Gene",
  remove_variants = "rs34637584"
)
fullSS_path <- example_fullSS()

Nalls23andMe_2019.results <- finemap_loci(
  top_SNPs = top_SNPs,
  results_dir = file.path(getwd(), "results"),
  loci = "BST1",
  dataset_name = "Nalls23andMe_2019",
  remove_tmps = FALSE,
  fullSS_path = fullSS_path,
  query_by = "tabix",
  snp_col = "RSID",
  pval_col = "p",
  effect_col = "beta",
  stderr_col = "se",
  freq_col = "freq",
  MAF_col = "calculate",
  bp_distance = 10000,
  min_MAF = 0.001,
  finemap_methods = c("FINEMAP", "SUSIE"),
  LD_reference = "UKB",
  download_method = "axel",
  plot.types = c()
)

Nalls23andMe_2019.results[SUSIE.CS > 0, list(SNP, FINEMAP.CS, FINEMAP.PP, SUSIE.CS, SUSIE.PP)]
readLines("results/GWAS/Nalls23andMe_2019/BST1/FINEMAP/data.snp", n = 3)
readLines("results/GWAS/Nalls23andMe_2019/BST1/FINEMAP/data.cred5")
