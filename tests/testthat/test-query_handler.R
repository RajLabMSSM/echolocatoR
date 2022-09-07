test_that("query_handler works", {

  fullSS_path <- echodata::example_fullSS()
  topSNPs <- echodata::topSNPs_Nalls2019
  dataset_type <- "GWAS"
  dataset_name <- "Nalls2019"
  locus <- "BST1"
  locus_dir <- echolocatoR:::construct_locus_dir(dataset_type = dataset_type,
                                                 dataset_name = dataset_name,
                                                 locus = locus)
  subset_path <- echolocatoR::: construct_subset_path(
    dataset_type = dataset_type,
    dataset_name = dataset_name,
    locus = locus)

  query <- echolocatoR::: query_handler(fullSS_path=fullSS_path,
                         colmap=echodata::construct_colmap(munged = TRUE),
                         locus_dir=locus_dir,
                         topSNPs=topSNPs,
                         subset_path=subset_path)
  main_cols <- c('SNP','CHR','BP','A1','A2',
  'FREQ','BETA','SE','P','N_CAS','N_CON')
  testthat::expect_true(methods::is(query,"data.table"))
  testthat::expect_equal(nrow(query), 2882)
  testthat::expect_true(all(main_cols %in% colnames(query)))
})
