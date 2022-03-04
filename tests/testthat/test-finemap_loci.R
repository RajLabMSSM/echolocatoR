test_that("finemap_loci works", {

  top_SNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")

  res <- echolocatoR::finemap_loci(
    fullSS_path = fullSS_path,
    top_SNPs = top_SNPs,
    loci = c("BST1","MEX3C"),
    dataset_name = "Nalls23andMe_2019",
    fullSS_genome_build = "hg19",
    zoom = c("1x","4x"),
    bp_distance = 250000,
    munged = TRUE)

  res
})
