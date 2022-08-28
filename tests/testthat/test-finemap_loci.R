test_that("finemap_loci works", {

  topSNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
  loci <- c("BST1","MEX3C")

  res <- echolocatoR::finemap_loci(
    fullSS_path = fullSS_path,
    topSNPs = topSNPs,
    loci = loci,
    dataset_name = "Nalls23andMe_2019",
    fullSS_genome_build = "hg19",
    finemap_methods = c("ABF","FINEMAP","SUSIE"),
    zoom = c("1x","4x"),
    bp_distance = 250000,
    munged = TRUE)

  testthat::expect_equal(names(res),c(loci,"merged_dat"))
  testthat::expect_gte(nrow(res$merged_dat), 2300)
  for(l in loci){
    testthat::expect_true(
      all(c("finemap_dat","locus_plot","LD_matrix",
            "LD_plot","locus_dir","arguments") %in% names(res[[l]]))
    )
  }
})
