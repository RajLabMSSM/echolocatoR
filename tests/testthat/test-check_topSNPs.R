test_that("multiplication works", {

  topSNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")

  #### From dataframe ####
  topSNPs2 <- check_topSNPs(topSNPs = topSNPs)
  testthat::expect_equal(topSNPs2$Locus, topSNPs$Locus)

  #### From fullSS ####
  fullSS <- data.table::fread(fullSS_path)
  fullSS$Locus <- paste0("Locus_",fullSS$CHR)
  data.table::fwrite(fullSS,fullSS_path)
  topSNPs3 <- check_topSNPs(topSNPs = "auto",
                            fullSS_path=fullSS_path)
  testthat::expect_equal(nrow(topSNPs3), length(unique(fullSS$Locus)))
})
