test_that("extract_snp_subset works", {

  subset_path <- echolocatoR:::construct_subset_path(dataset_type = "GWAS",
                                                     dataset_name = "Nalls",
                                                     locus = "BST1")
  fullSS_path <- echodata::example_fullSS()
  topSNPs <- echodata::topSNPs_Nalls2019

  query <- echolocatoR:::extract_snp_subset(
    subset_path=subset_path,
    colmap=echodata::construct_colmap(munged = TRUE),
    fullSS_path=fullSS_path,
    topSNPs=topSNPs,
    LD_reference="1KGphase3")

  testthat::expect_true(methods::is(query,"data.table"))
  testthat::expect_equal(nrow(query), 2882)


  #### Restore save query ####
  query <- echolocatoR:::extract_snp_subset(
    subset_path=subset_path,
    colmap=echodata::construct_colmap(munged = TRUE),
    fullSS_path=fullSS_path,
    topSNPs=topSNPs,
    LD_reference="1KGphase3")
  testthat::expect_true(methods::is(query,"data.table"))
  testthat::expect_equal(nrow(query), 2882)

})
