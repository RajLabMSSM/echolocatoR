test_that("check_deprecated works", {

  topSNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS()

  ## Test that deprecated colmap args (chrom_col, position_col) raise an error
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

  ## Test that deprecated non-colmap args (PAINTOR_QTL_datasets, server)

  ## raise a warning via check_deprecated directly,
  ## without running the full pipeline (which requires gcc/FINEMAP).
  fake_call <- quote(finemap_loci(
    fullSS_path = "test",
    topSNPs = "test",
    loci = "BST1",
    PAINTOR_QTL_datasets = c("test"),
    server = TRUE
  ))
  testthat::expect_warning(
    echolocatoR::check_deprecated(
      fun = "finemap_loci",
      args = fake_call
    )
  )

  fun <- function(){
    match.call()
  }
})
