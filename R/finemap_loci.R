#' Fine-map multiple loci
#'
#' \pkg{echolocatoR} will automatically fine-map each locus.
#' Uses the \code{topSNPs} data.frame to define locus coordinates.
#'
#' @family MAIN
#' @param loci The list of loci you want to fine-map.
#' If \code{subset_path="auto"} (\emph{default}),
#'  a locus subset file name is automatically constructed as:
#' \emph{Data/<dataset_type>/<dataset_name>/<locus>/
#' Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}
#' @param loci Character list of loci in \strong{Locus} col of \code{topSNPs}.
#' @param return_all Return a nested list of various the pipeline's outputs
#' including plots, tables, and file paths (default: \code{TRUE}).
#' If \code{FALSE}, instead only returns a single merged
#'  \link[data.table]{data.table} containing the results from all loci.
#' @param use_tryCatch If an error is encountered in one locus,
#' the pipeline will continue to try running the rest of the loci
#'  (default: \code{use_tryCatch=TRUE}). This avoid stopping all analyses due
#'  to errors that only affect some loci,
#'  but currently prevents debugging via traceback.
#' @inheritParams finemap_locus
#' @inheritParams echoconda::activate_env
#' @inheritParams echodata::filter_snps
#' @inheritParams echoLD::get_LD
#' @inheritParams echoLD::filter_LD
#' @inheritParams echoplot::plot_locus
#' @inheritParams echofinemap::multifinemap
#' @inheritParams echodata::standardize
#' @inheritParams echodata::find_consensus_snps
#'
#' @return A merged data.frame with all fine-mapping results from all loci.
#'
#' @export
#' @importFrom echodata find_consensus_snps construct_colmap gene_locus_list
#' @importFrom data.table rbindlist data.table
#' @importFrom echofinemap POLYFUN_install
#'
#' @examples
#' topSNPs <- echodata::topSNPs_Nalls2019
#' fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
#'
#' res <- echolocatoR::finemap_loci(
#'   fullSS_path = fullSS_path,
#'   topSNPs = topSNPs,
#'   loci = c("BST1","MEX3C"),
#'   finemap_methods = c("ABF","FINEMAP","SUSIE"),
#'   dataset_name = "Nalls23andMe_2019",
#'   fullSS_genome_build = "hg19",
#'   bp_distance = 10000,
#'   munged = TRUE)
finemap_loci <- function(#### Main args ####
                         loci = NULL,
                         fullSS_path,
                         fullSS_genome_build = NULL,
                         results_dir = file.path(tempdir(),"results"),
                         dataset_name = "dataset_name",
                         dataset_type = "GWAS",
                         topSNPs = "auto",
                         #### Force new args ####
                         force_new_subset = FALSE,
                         force_new_LD = FALSE,
                         force_new_finemap = FALSE,
                         #### Fine-mapping args ####
                         finemap_methods = c("ABF","FINEMAP","SUSIE"),
                         finemap_args = NULL,
                         n_causal = 5,
                         credset_thresh = .95,
                         consensus_thresh = 2,
                         fillNA = 0,
                         conditioned_snps = "auto",
                         priors_col = NULL,
                         #### Colname mapping args ####
                         munged = FALSE,
                         colmap = echodata::construct_colmap(munged = munged),
                         compute_n = "ldsc",
                         #### LD args ####
                         LD_reference = "1KGphase3",
                         LD_genome_build = "hg19",
                         leadSNP_LD_block = FALSE,
                         superpopulation = "EUR",
                         download_method = "axel",
                         #### SNP filter args ####
                         bp_distance = 500000,
                         min_POS = NA,
                         max_POS = NA,
                         min_MAF = NA,
                         trim_gene_limits = FALSE,
                         max_snps = NULL,
                         min_r2 = 0,
                         remove_variants = FALSE,
                         remove_correlates = FALSE,
                         #### Misc args ####
                         query_by = "tabix",
                         case_control = TRUE,
                         qtl_prefixes = NULL,
                         #### PLotting args ####
                         plot_types = c("simple"),
                         zoom = "1x",
                         nott_epigenome = FALSE,
                         nott_show_placseq = FALSE,
                         nott_binwidth = 200,
                         nott_bigwig_dir = NULL,
                         xgr_libnames = NULL,
                         roadmap = FALSE,
                         roadmap_query = NULL,
                         #### General args ####
                         remove_tmps = TRUE,
                         conda_env = "echoR_mini",
                         return_all = TRUE,
                         use_tryCatch = TRUE,
                         seed = 2022,
                         nThread = 1,
                         verbose = TRUE,
                         #### Deprecated args ####
                         top_SNPs = deprecated(),
                         PP_threshold = deprecated(),
                         consensus_threshold = deprecated(),
                         plot.Nott_epigenome = deprecated(),
                         plot.Nott_show_placseq = deprecated(),
                         plot.Nott_binwidth = deprecated(),
                         plot.Nott_bigwig_dir = deprecated(),
                         plot.Roadmap = deprecated(),
                         plot.Roadmap_query = deprecated(),
                         plot.XGR_libnames = deprecated(),
                         server = deprecated(),
                         plot.types = deprecated(),
                         plot.zoom = deprecated(),
                         QTL_prefixes = deprecated(),
                         vcf_folder = deprecated(),
                         probe_path = deprecated(),
                         file_sep=deprecated(),
                         ## Deprecated colmap args
                         chrom_col=deprecated(),
                         chrom_type=deprecated(),
                         position_col=deprecated(),
                         snp_col=deprecated(),
                         pval_col=deprecated(),
                         effect_col=deprecated(),
                         stderr_col=deprecated(),
                         tstat_col=deprecated(),
                         locus_col=deprecated(),
                         freq_col=deprecated(),
                         MAF_col=deprecated(),
                         A1_col = deprecated(),
                         A2_col = deprecated(),
                         gene_col=deprecated(),
                         N_cases_col=deprecated(),
                         N_controls_col=deprecated(),
                         N_cases=deprecated(),
                         N_controls=deprecated(),
                         proportion_cases=deprecated(),
                         sample_size=deprecated(),
                         PAINTOR_QTL_datasets=deprecated()
                         ){
  # echoverseTemplate:::source_all();
  # echoverseTemplate:::args2vars(finemap_loci);

  check_deprecated(fun="finemap_loci",
                   args=match.call(call = sys.call(sys.parent(2))))
  #### Check PolyFun/PAINTOR are installed early on (if needed) ####
  if(any(grepl("polyfun",finemap_methods, ignore.case = TRUE))){
    echofinemap::POLYFUN_install()
  }
  if(any(grepl("PAINTOR",finemap_methods, ignore.case = TRUE))){
    echofinemap:::PAINTOR_install()
  }
  #### Parallise data.table functions ####
  data.table::setDTthreads(threads = nThread);
  #### Get loci (if not supplied) ####
  loci <- echodata::gene_locus_list(loci = loci,
                                    topSNPs = topSNPs,
                                    dataset_type = dataset_type,
                                    verbose = verbose)
  #### Validate fullSS genome build ####
  fullSS_genome_build <- check_genome(fullSS_genome_build = fullSS_genome_build,
                                      munged = colmap$munged,
                                      fullSS_path = fullSS_path)
  #### Select SNPs to condition on (for COJO) ####
  conditioned_snps <- snps_to_condition(conditioned_snps = conditioned_snps,
                                        topSNPs = topSNPs,
                                        loci = loci)

  #### Iterate fine-mapping over loci ####
  t1 <- Sys.time()
  FINEMAP_DAT <- lapply(seq_len(length(loci)),
                          function(i){
    t1_locus <- Sys.time()
    finemap_dat <- NULL
    locus <- loci[i]
    #### Locus header ####
    cli::cat_boxx(
      cli::col_br_cyan(
        paste(cli::col_br_white(cli::style_blurred(")))>")),
              bat_icon(), locus,
              paste0("[locus ",i," / ",length(loci),"]"),
              bat_icon(),
              cli::col_br_white(cli::style_blurred("<((("))
              )
      )
    )
    tryCatch({
      gene_limits <- arg_list_handler(trim_gene_limits, i)
      conditioned_snp <- arg_list_handler(conditioned_snps, i)
      min_pos <- arg_list_handler(min_POS, i)
      max_pos <- arg_list_handler(max_POS, i)
      LD_ref <- arg_list_handler(LD_reference, i)

      out_list <- finemap_locus(locus=locus,
                                topSNPs=topSNPs,
                                fullSS_path=fullSS_path,
                                fullSS_genome_build=fullSS_genome_build,
                                munged=munged,
                                colmap=colmap,
                                compute_n=compute_n,
                                priors_col=priors_col,

                                LD_genome_build=LD_genome_build,
                                results_dir=results_dir,
                                finemap_methods=finemap_methods,
                                finemap_args=finemap_args,
                                force_new_subset=force_new_subset,
                                force_new_LD=force_new_LD,
                                force_new_finemap=force_new_finemap,
                                dataset_name=dataset_name,
                                dataset_type=dataset_type,
                                n_causal=n_causal,
                                bp_distance=bp_distance,

                                 LD_reference=LD_ref,
                                 superpopulation=superpopulation,
                                 download_method=download_method,
                                 min_POS=min_pos,
                                 max_POS=max_pos,
                                 min_MAF=min_MAF,

                                 trim_gene_limits=gene_limits,
                                 max_snps=max_snps,
                                 min_r2=min_r2,
                                 leadSNP_LD_block=leadSNP_LD_block,
                                 query_by=query_by,
                                 remove_variants=remove_variants,
                                 remove_correlates=remove_correlates,
                                 conditioned_snps=conditioned_snps,
                                 remove_tmps=remove_tmps,
                                 plot.types=plot.types,
                                 credset_thresh=credset_thresh,
                                 consensus_thresh=consensus_thresh,
                                 case_control=case_control,
                                 qtl_prefixes=qtl_prefixes,

                                 plot_types=plot_types,
                                 zoom=zoom,
                                 nott_epigenome=nott_epigenome,
                                 nott_show_placseq=nott_show_placseq,
                                 nott_binwidth=nott_binwidth,
                                 nott_bigwig_dir=nott_bigwig_dir,
                                 xgr_libnames=xgr_libnames,
                                 roadmap=roadmap,
                                 roadmap_query=roadmap_query,

                                 seed=seed,
                                 conda_env=conda_env,
                                 nThread=nThread,
                                 verbose=verbose)
      #### Output list reminder ####
      # list(finemap_dat
      #     locus_plot
      #     LD_matrix
      #     LD_plot
      #     locus_dir
      #     arguments)
      messager("Formatting locus results.",v=verbose)
      if(return_all) return(out_list)
      finemap_dat <- out_list$finemap_dat
      if(!"Locus" %in% colnames(finemap_dat)){
        finemap_dat <- data.table::data.table(Locus=locus,
                                              finemap_dat)
      }
      cat('  \n')
    }, error = set_tryCatch(use_tryCatch = use_tryCatch)) ## end tryCatch()
    report_time(t1 = t1_locus,
                prefix = paste("Locus",locus,"complete in:"))
    return(finemap_dat)
  }) # end loop
  #### Prepare results to return ####
  steps("postprocess")
  FINEMAP_DAT <- postprocess_data(FINEMAP_DAT=FINEMAP_DAT,
                                  loci=loci,
                                  credset_thresh=credset_thresh,
                                  consensus_thresh=consensus_thresh,
                                  return_all=return_all,
                                  verbose=verbose)
  report_time(t1 = t1,
              prefix = "All loci done in:")
  return(FINEMAP_DAT)
}
