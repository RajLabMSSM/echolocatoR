test_that("check_genome works", {

  fullSS_path <- echodata::example_fullSS()
  #### Infer ###
  ## Requires several large databases to be installed first. Skip.
  # build <- check_genome(fullSS_path = fullSS_path,
  #                       fullSS_genome_build = "hg19",
  #                       sampled_snps = 100,
  #                       munged = TRUE)
  #### Standardize and pass back ####
  build1 <- check_genome(fullSS_genome_build="hg19")
  testthat::expect_equal(build1,"GRCH37")
  build2 <- check_genome(fullSS_genome_build="grcH38")
  testthat::expect_equal(build2,"GRCH38")
})
