test_that("finemap_loci works", {

  topSNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
  loci <- c("BST1","MEX3C")

  run_tests <- function(res){
    for(l in loci){
      testthat::expect_true(
        all(c("finemap_dat","locus_plot","LD_matrix",
              "LD_plot","locus_dir","arguments") %in% names(res[[l]]))
      )
    }
  }

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
  run_tests(res = res)

  #### When use_tryCatch=TRUE ####
  res <- echolocatoR::finemap_loci(
    fullSS_path = fullSS_path,
    topSNPs = topSNPs,
    loci = c("typoooo","MEX3C"),
    finemap_methods = c("ABF","FINEMAP","SUSIE"),
    dataset_name = "Nalls23andMe_2019",
    fullSS_genome_build = "hg19",
    bp_distance = 10000,
    use_tryCatch = TRUE,
    munged = TRUE)
  testthat::expect_null(res$typoooo)
  testthat::expect_true(methods::is(res$MEX3C$finemap_dat,"data.table"))

  #### When use_tryCatch=FALSE ####
  testthat::expect_error(
    res <- echolocatoR::finemap_loci(
      fullSS_path = fullSS_path,
      topSNPs = topSNPs,
      loci = c("typoooo","MEX3C"),
      finemap_methods = c("ABF","FINEMAP","SUSIE"),
      dataset_name = "Nalls23andMe_2019",
      fullSS_genome_build = "hg19",
      bp_distance = 10000,
      use_tryCatch = FALSE,
      munged = TRUE)
  )


  #### Using custom LD: rds + vcf ####
  ## BST1
  ld1 <- system.file("extdata", "BST1.1KGphase3.vcf.bgz",
      package = "echodata")
  ## MEX3C (fake)
  ld2 <- tempfile(fileext = ".tsv.gz")
  ld_mat <- echodata::BST1_LD_matrix
  mex3c <- echodata::MEX3C[seq_len(nrow(ld_mat)),]
  ld_mat <- ld_mat |> `rownames<-`(mex3c$SNP) |> `colnames<-`(mex3c$SNP)
  data.table::fwrite(data.table::data.table(ld_mat,keep.rownames = "rsid"),
                     ld2,
                     sep="\t")

  LD_reference <- list(ld1, ld2)
  res <- echolocatoR::finemap_loci(
    fullSS_path = fullSS_path,
    topSNPs = topSNPs,
    loci = c("BST1","MEX3C"),
    finemap_methods = c("ABF","FINEMAP","SUSIE"),
    dataset_name = "Nalls23andMe_2019",
    fullSS_genome_build = "hg19",
    LD_reference = LD_reference,
    bp_distance = 10000,
    munged = TRUE)
  run_tests(res = res)
})
