test_that("finemap_locus works", {

  topSNPs <- echodata::topSNPs_Nalls2019
  fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")

  res <- echolocatoR::finemap_locus(
    fullSS_path = fullSS_path,
    topSNPs = topSNPs,
    # results_dir = "/Desktop/res",
    locus = "BST1",
    dataset_name = "Nalls2019",
    fullSS_genome_build = "hg19",
    zoom = c("1x","4x"),
    bp_distance = 25000,
    n_causal = 5,
    force_new_finemap = TRUE,
    plot_types = c("simple","fancy","LD"),
    roadmap = TRUE,
    roadmap_query = "E053",
    # nott_epigenome = TRUE,
    # nott_show_placseq = TRUE,
    munged = TRUE)

  testthat::expect_true(
    all(c("finemap_dat","locus_plot","LD_matrix",
          "LD_plot","locus_dir","arguments") %in% names(res))
  )

  testthat::expect_equal(nrow(res$finemap_dat), 210)
  testthat::expect_equal(nrow(res$LD_matrix), 210)
  # failing atm (due to roadmap error in fancy plot)
  testthat::expect_equal(names(res$locus_plot),c("simple","fancy"))
  testthat::expect_equal(names(res$locus_plot$simple),c("1x","4x"))
  testthat::expect_gte(
      sum(echodata::BST1$SNP %in% res$locus_plot$simple$`1x`$data$SNP),
      200
  )
  testthat::expect_true(
    all(c("FINEMAP","LD","Multi-finemap",
          "multiview.BST1.1KGphase3.1x.png") %in% list.files(res$locus_dir))
  )
})
