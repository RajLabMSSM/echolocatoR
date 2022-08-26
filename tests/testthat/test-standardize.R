test_that("standardize_subset works", {

  BST1 <- echodata::BST1
  #### Screw up Freq to see if function can fix it and infer MAF ####
  BST1$rsid <- BST1$SNP
  BST1 <- data.frame(BST1)[,!colnames(BST1) %in% c("MAF","SNP")]
  BST1[c(10,30,55),"Freq"] <- 0
  BST1[c(12,22),"Freq"] <- NA

  subset_path <- file.path(tempdir(),"BST1.tsv")
  data.table::fwrite(BST1, subset_path)
  query_mod <- echolocatoR:::standardize_subset(subset_path=subset_path,
                                                locus="BST1",
                                                snp_col = "rsid",
                                                MAF_col="calculate")
  testthat::expect_equal(
    colnames(query_mod),
    c('CHR','POS','SNP','P','Effect','StdErr','A1',
      'A2','Freq','MAF','N_cases','N_controls',
      'proportion_cases','N','t_stat','leadSNP')
  )
  testthat::expect_equal(dim(query_mod),c(6216,16))

})
