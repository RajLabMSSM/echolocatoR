test_that("check_deprecated works", {

  topSNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS()

  testthat::expect_error(
    echolocatoR::finemap_loci(
      fullSS_path = fullSS_path,
      topSNPs = topSNPs,
      loci = c("BST1","MEX3C"),
      munged = TRUE,
      fullSS_genome_build = "hg19",
      bp_distance=10000,
      chrom_col = "CHR",
      position_col = "BP")
  )
  testthat::expect_warning(
    echolocatoR::finemap_loci(
      fullSS_path = fullSS_path,
      topSNPs = topSNPs,
      loci = c("BST1"),
      munged = TRUE,
      fullSS_genome_build = "hg19",
      bp_distance=10000,
      PAINTOR_QTL_datasets = c("test"),
      server = TRUE
      )
  )
})
