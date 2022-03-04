#' Fine-map multiple loci
#'
#' \pkg{echolocatoR} will automatically fine-map each locus.
#' Uses the \code{top_SNPs} data.frame to define locus coordinates.
#'
#' @family MAIN
#' @param loci The list of loci you want to fine-map.
#' If \code{subset_path="auto"} (\emph{default}),
#'  a locus subset file name is automatically constructed as:
#' \emph{Data/<dataset_type>/<dataset_name>/<locus>/
#' Multi-finemap/<locus>_<dataset_name>_Multi-finemap.tsv.gz}
#' @inheritParams finemap_locus
#' @inheritParams echoconda::activate_env
#' @inheritParams echodata::filter_snps
#' @inheritParams echoLD::load_or_create
#' @inheritParams echoLD::filter_LD
#' @inheritParams echoplot::plot_locus
#'
#' @return A merged data.frame with all fine-mapping results from all loci.
#'
#' @export
#' @importFrom echodata MUNGESUMSTATS.col_map find_consensus_snps
#' @importFrom stats setNames
#' @importFrom data.table rbindlist data.table
#'
#' @examples
#' top_SNPs <- echodata::topSNPs_Nalls2019
#' fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
#'
#' Nalls23andMe_2019.results <- echolocatoR::finemap_loci(
#'   fullSS_path = fullSS_path,
#'   top_SNPs = top_SNPs,
#'   loci = c("BST1","MEX3C"),
#'   dataset_name = "Nalls23andMe_2019",
#'   fullSS_genome_build = "hg19",
#'   bp_distance = 250000,
#'   munged = TRUE)
finemap_loci <- function(#### Main args ####
                         loci,
                         fullSS_path,
                         fullSS_genome_build=NULL,
                         results_dir=file.path(tempdir(),"results"),
                         dataset_name="dataset_name",
                         dataset_type="GWAS",
                         top_SNPs="auto",
                         #### Force new args ####
                         force_new_subset=FALSE,
                         force_new_LD=FALSE,
                         force_new_finemap=TRUE,
                         #### Fine-mapping args ####
                         finemap_methods=c("ABF","FINEMAP",
                                           "SUSIE","POLYFUN_SUSIE"),
                         finemap_args=NULL,
                         n_causal=5,
                         PP_threshold=.95,
                         consensus_threshold=2,
                         fillNA=0,
                         conditioned_snps = "auto",
                         #### Colname mapping args ####
                         munged = FALSE,
                         chrom_col="CHR",
                         chrom_type=NULL,
                         position_col="POS",
                         snp_col="SNP",
                         pval_col="P",
                         effect_col="Effect",
                         stderr_col="StdErr",
                         tstat_col="t-stat",
                         locus_col="Locus",
                         freq_col="Freq",
                         MAF_col="MAF",
                         A1_col = "A1",
                         A2_col = "A2",
                         gene_col="Gene",
                         N_cases_col="N_cases",
                         N_controls_col="N_controls",
                         N_cases=NULL,
                         N_controls=NULL,
                         proportion_cases="calculate",
                         sample_size=NULL,
                         #### LD args ####
                         LD_reference="1KGphase1",
                         LD_genome_build="hg19",
                         leadSNP_LD_block=FALSE,
                         superpopulation="EUR",
                         download_method="axel",
                         #### SNP filter args ####
                         bp_distance=500000,
                         min_POS=NA,
                         max_POS=NA,
                         min_MAF=NA,
                         trim_gene_limits=FALSE,
                         max_snps=NULL,
                         min_r2=0,
                         remove_variants=FALSE,
                         remove_correlates=FALSE,
                         #### Misc args ####
                         query_by="tabix",
                         PAINTOR_QTL_datasets=NULL,
                         case_control=TRUE,
                         qtl_prefixes=NULL,
                         #### PLotting args ####
                         plot_types=c("simple"),
                         zoom="1x",
                         nott_epigenome=FALSE,
                         nott_show_placseq=FALSE,
                         nott_binwidth=200,
                         nott_bigwig_dir=NULL,
                         xgr_libnames=NULL,
                         roadmap=FALSE,
                         roadmap_query=NULL,
                         #### General args ####
                         remove_tmps=TRUE,
                         conda_env="echoR",
                         return_all=TRUE,
                         nThread=1,
                         verbose=TRUE,
                         #### Deprecated args ####
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
                         file_sep=deprecated()){
  #### Conda env setup ####
  if(tolower(conda_env)=="echor"){
    conda_env <- echoconda::yaml_to_env(yaml_path = conda_env,
                                        verbose = FALSE)
  }
  echoconda::activate_env(conda_env = conda_env,
                          verbose = verbose)
  #### Check if CLI tools available ####
  ## These CLI tools are no longer essential as I've replaced
  ## them all with R-native alternatives.
  check_tools(check_tabix=TRUE,
              stop_for_tabix=FALSE,
              check_bcftools=TRUE,
              #stop_for_bcftools = LD_reference %in% c("1KGphase1","1KGphase3"),
              stop_for_bcftools=FALSE,
              conda_env=conda_env)
  fullSS_genome_build <- check_genome(gbuild=fullSS_genome_build,
                                      munged=munged,
                                      fullSS_path=fullSS_path)
  data.table::setDTthreads(threads = nThread);
  conditioned_snps <- snps_to_condition(conditioned_snps, top_SNPs, loci);

  if(detect_genes(loci = loci, verbose  = FALSE)){
    messager("Reassigning gene-specific locus names",v=verbose)
    loci <- stats::setNames(paste(unname(loci),names(loci),sep="_"),
                            names(loci))
  }
  #### Iterate over loci ####
  FINEMAP_DAT <- lapply(seq_len(length(unique(loci))), function(i){
    start_gene <- Sys.time()
    finemap_dat <- NULL
    locus <- loci[i]
    messager("\n)   )  ) ))))))}}}}}}}}",
             locus," (",i ," / ",length(loci),")",
             "{{{{{{{{{(((((( (  (   (");
    try({
      # lead_SNP <- arg_list_handler(conditioned_snps, i)
      gene_limits <- arg_list_handler(trim_gene_limits, i)
      conditioned_snp <- arg_list_handler(conditioned_snps, i)
      min_pos <- arg_list_handler(min_POS, i)
      max_pos <- arg_list_handler(max_POS, i)
      LD_ref <- arg_list_handler(LD_reference, i)

      out_list <- finemap_locus(locus=locus,
                                   top_SNPs=top_SNPs,
                                   fullSS_path=fullSS_path,
                                   fullSS_genome_build=fullSS_genome_build,
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

                                   chrom_col=chrom_col,
                                   chrom_type=chrom_type,
                                   position_col=position_col,
                                   snp_col=snp_col,
                                   pval_col=pval_col,
                                   effect_col=effect_col,
                                   stderr_col=stderr_col,
                                   tstat_col=tstat_col,
                                   locus_col=locus_col,
                                   MAF_col=MAF_col,
                                   freq_col=freq_col,
                                   A1_col=A1_col,
                                   A2_col=A2_col,
                                   gene_col=gene_col,
                                   N_cases_col=N_cases_col,
                                   N_controls_col=N_controls_col,
                                   N_cases=N_cases,
                                   N_controls=N_controls,
                                   proportion_cases=proportion_cases,
                                   sample_size=sample_size,

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
                                   PAINTOR_QTL_datasets=PAINTOR_QTL_datasets,
                                   PP_threshold=PP_threshold,
                                   consensus_threshold=consensus_threshold,
                                   case_control=case_control,
                                   qtl_prefixes=qtl_prefixes,

                                   zoom=zoom,
                                   nott_epigenome=nott_epigenome,
                                   nott_show_placseq=nott_show_placseq,
                                   nott_binwidth=nott_binwidth,
                                   nott_bigwig_dir=nott_bigwig_dir,
                                   xgr_libnames=xgr_libnames,
                                   roadmap=roadmap,
                                   roadmap_query=roadmap_query,

                                   munged=munged,
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
      if(return_all) return(out_list)
      finemap_dat <- out_list$finemap_dat
      if(!"Locus" %in% colnames(finemap_dat)){
        finemap_dat <- data.table::data.table(Locus=locus,
                                              finemap_dat)
      }
      cat('  \n')
    }) ## end try()
    end_gene <- Sys.time()
    message("Fine-mapping complete in:")
    print(round(end_gene-start_gene,1))
    return(finemap_dat)
  }) # end for loop
  #### Prepare results to return ####
  FINEMAP_DAT <- postprocess_data(FINEMAP_DAT=FINEMAP_DAT,
                                  loci=loci,
                                  return_all=return_all,
                                  verbose=verbose)
  return(FINEMAP_DAT)
}
