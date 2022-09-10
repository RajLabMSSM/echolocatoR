#' Run \pkg{echolocatoR} pipeline on a single locus
#'
#' Unlike \code{finemap_loci}, you don't need to provide a \code{topSNPs}
#' data.frame. Instead, just manually provide the coordinates of the locus
#' you want to fine-map.
#'
#' The primary functions of \pkg{echolocatoR} that expedite fine-mapping
#'  by wrapping many other \pkg{echolocatoR} functions into one.
#'  Encompasses steps including:
#'  \describe{
#'  \item{Subset & standardize}{Extract subsets of the full summary stats
#'   GWAS or QTL file and reformat them to be compatible with
#'   \pkg{echolocatoR}'s various functions }
#'  \item{Calculate linkage disequilibrium}{Download and prepare the
#'  necessary LD matrix.}
#'  \item{Fine-map}{Run various fine-mapping tools and merge the results
#'   into a single multi-finemap data.frame.}
#'  \item{Plot}{Summarise the results in a multi-track plot for each locus.}
#'  }
#'
#' @param locus Locus name to fine-map (e.g. \code{"BIN1"}).
#' Can be named to indicate a specific gene within a QTL locus
#' (e.g. \code{c(ENSG00000136731="BIN1")}).
#' @param fullSS_path Path to the full summary statistics file (GWAS or QTL)
#' that you want to fine-map.
#' It is usually best to provide the absolute path rather
#' than the relative path.
#' @param fullSS_genome_build Genome build of the full summary statistics
#'  (\code{fullSS_path}). Can be "GRCH37" or "GRCH38" or one of their synonyms..
#' If \code{fullSS_genome_build==NULL} and \code{munged=TRUE},
#' infers genome build (hg19 vs. hg38)
#' from summary statistics using \link[MungeSumstats]{get_genome_builds}.
#' @param LD_genome_build Genome build of the LD panel.
#' This is automatically assigned to the correct genome build for each
#' LD panel except when the user supplies custom vcf/LD files.
#' @param results_dir Where to store all results.
#' \strong{IMPORTANT!:} It is usually best to provide the absolute path
#' rather than the relative path.
#' This is especially important for \emph{FINEMAP}.
#' @param query_by Choose which method you want to use to extract
#'  locus subsets from the full summary stats file.
#' Methods include:
#' \describe{
#' \item{"tabix"}{Convert the full summary stats file in an indexed tabix file.
#'  Makes querying lightning fast after the initial conversion is done.
#'   (\emph{default})}
#' \item{"coordinates"}{Extract locus subsets using min/max genomic
#' coordinates with \emph{awk}.}
#' }
#' @param dataset_name The name you want to assign to the dataset
#' being fine-mapped,
#' This will be used to name the subdirectory where your
#' results will be stored
#' (e.g. \emph{Data/GWAS/<dataset_name>}).
#' Don't use special characters (e.g.".", "/").
#' @param dataset_type The kind dataset you're fine-mapping
#' (e.g. GWAS, eQTL, tQTL).
#' This will also be used when creating the subdirectory where your results
#' will be stored
#' (e.g. \emph{Data/<dataset_type>/Kunkle_2019}).
#' @param topSNPs A data.frame with the genomic coordinates of the lead SNP
#' for each locus.
#' The lead SNP will be used as the center of the window when extracting
#' subset from the full GWAS/QTL summary statistics file.
#' Only one SNP per \strong{Locus} should be included.
#' At minimum, \code{topSNPs} should include the following columns:
#' \describe{
#' \item{\emph{Locus}}{A unique name for each locus. Often,
#'  loci are named after a relevant gene (e.g. LRRK2) or based on
#'   the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
#' \item{\emph{CHR}}{The chromosome that the SNP is on.
#'  Can be "chr12" or "12" format.}
#' \item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
#' }
#'
#' @param force_new_subset By default, if a subset of the full
#'  summary stats file for a given locus is already present,
#' then \pkg{echolocatoR} will just use the pre-existing file.
#' Set \code{force_new_subset=T} to override this and extract a new subset.
#' Subsets are saved in the following path structure:
#' \emph{Data/\<dataset_type\>/\<dataset_name\>/\<locus\>/Multi-finemap/
#' \<locus\>_\<dataset_name\>_Multi-finemap.tsv.gz}
#' @param force_new_finemap By default, if an fine-mapping results file for
#'  a given locus is already present,
#' then \pkg{echolocatoR} will just use the preexisting file.
#' Set \code{force_new_finemap=T} to override this and re-run fine-mapping.
#'
#' @param finemap_methods Which fine-mapping methods you want to use.
#' @param n_causal The maximum number of potential causal SNPs per locus.
#' This parameter is used somewhat differently by different fine-mapping tools.
#' See tool-specific functions for details.
#'
#' @param munged Whether \code{fullSS_path} have already been
#' standardised/filtered full summary stats
#' with \link[MungeSumstats]{format_sumstats}.
#' If \code{munged=FALSE} you'll need to provide the necessary
#'  column names to the \code{colmap} argument.
#' @param colmap Column name mappings in in \code{fullSS_path}. Must be a named
#' list. Can use \link[echodata]{construct_colmap} to assist with this. This
#' function can be used in two different ways:
#' \itemize{
#' \item{\code{munged=FALSE} : }{When \code{munged=FALSE},
#'  you will need to provide the necessary column names to the
#'  \code{colmap} argument (\emph{default}).}
#'  \item{\code{munged=TRUE} : }{ Alternatively, instead of filling out
#'  each argument in
#' \link[echodata]{construct_colmap}, you can simply set \code{munged=TRUE}
#'  if  \code{fullSS_path} has already been munged with
#'  \link[MungeSumstats]{format_sumstats}.
#'  }
#' }
#' @param conditioned_snps Which SNPs to conditions on when fine-mapping
#' with \emph{COJO}.
#' @param plot_types Which kinds of plots to include.
#' Options:
#' \itemize{
#' \item{"simple"}{Just plot the following tracks: GWAS,
#' fine-mapping, gene models}
#' \item{"fancy"}{Additionally plot XGR annotation tracks
#' (XGR, Roadmap, Nott2019).}
#' ' \item{"LD"}{LD heatmap showing the 10 SNPs surrounding the lead SNP.}
#' }
#' @param remove_tmps Whether to remove any temporary files
#'  (e.g. FINEMAP output files) after the pipeline is done running.
#' @param seed Set the seed for all functions where this is possible.
#'
#' @param case_control [deprecated]
#' @param top_SNPs [deprecated]
#' @param PP_threshold [deprecated]
#' @param top_SNPs [deprecated]
#' @param consensus_threshold [deprecated]
#' @param plot.Nott_epigenome [deprecated]
#' @param plot.Nott_show_placseq [deprecated]
#' @param plot.Nott_binwidth [deprecated]
#' @param plot.Nott_bigwig_dir [deprecated]
#' @param plot.XGR_libnames [deprecated]
#' @param plot.Roadmap [deprecated]
#' @param plot.Roadmap_query [deprecated]
#' @param server [deprecated]
#' @param plot.types [deprecated]
#' @param plot.zoom [deprecated]
#' @param QTL_prefixes [deprecated]
#' @param vcf_folder [deprecated]
#' @param probe_path [deprecated]
#' @param file_sep [deprecated]
#' @param chrom_col [deprecated]
#' @param chrom_type [deprecated]
#' @param position_col [deprecated]
#' @param freq_col [deprecated]
#' @param snp_col [deprecated]
#' @param pval_col [deprecated]
#' @param effect_col [deprecated]
#' @param stderr_col [deprecated]
#' @param tstat_col [deprecated]
#' @param locus_col [deprecated]
#' @param MAF_col [deprecated]
#' @param A1_col [deprecated]
#' @param A2_col [deprecated]
#' @param gene_col [deprecated]
#' @param N_cases_col [deprecated]
#' @param N_controls_col [deprecated]
#' @param N_cases [deprecated]
#' @param N_controls [deprecated]
#' @param proportion_cases [deprecated]
#' @param sample_size [deprecated]
#' @param PAINTOR_QTL_datasets [deprecated]
#'
#' @family MAIN
#' @inheritParams echoconda::activate_env
#' @inheritParams echodata::filter_snps
#' @inheritParams echoLD::get_LD
#' @inheritParams echoLD::filter_LD
#' @inheritParams echoplot::plot_locus
#' @inheritParams echofinemap::multifinemap
#' @inheritParams echodata::standardize
#' @inheritParams echodata::find_consensus_snps
#'
#' @importFrom echodata construct_colmap
#' @importFrom echoplot plot_locus
#' @importFrom echoLD get_LD filter_LD subset_common_snps
#' @importFrom echofinemap multifinemap
#' @export
#' @examples
#' topSNPs <- echodata::topSNPs_Nalls2019
#' fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
#'
#' res <- echolocatoR::finemap_locus(
#'   fullSS_path = fullSS_path,
#'   topSNPs = topSNPs,
#'   locus = "BST1",
#'   finemap_methods = c("ABF","FINEMAP","SUSIE"),
#'   dataset_name = "Nalls23andMe_2019",
#'   fullSS_genome_build = "hg19",
#'   bp_distance = 1000,
#'   munged = TRUE)
finemap_locus <- function(#### Main args ####
                          locus,
                          fullSS_path,
                          fullSS_genome_build = NULL,
                          results_dir = file.path(tempdir(),"results"),
                          dataset_name = "dataset_name",
                          dataset_type = "GWAS",
                          case_control = TRUE,
                          topSNPs = "auto",
                          #### Force new args ####
                          force_new_subset = FALSE,
                          force_new_LD = FALSE,
                          force_new_finemap = FALSE,
                          #### Fine-mapping args ####
                          finemap_methods = c("ABF","FINEMAP","SUSIE"),
                          finemap_args=NULL,
                          n_causal=5,
                          credset_thresh=.95,
                          consensus_thresh=2,
                          fillNA=0,
                          conditioned_snps=NULL,
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
                          qtl_suffixes = NULL,
                          #### PLotting args ####
                          plot_types = c("simple"),
                          zoom = "1x",
                          show_plot = TRUE,
                          tx_biotypes = NULL,
                          nott_epigenome = FALSE,
                          nott_show_placseq = FALSE,
                          nott_binwidth = 200,
                          nott_bigwig_dir = NULL,
                          xgr_libnames = NULL,
                          roadmap = FALSE,
                          roadmap_query = NULL,
                          #### General args ####
                          remove_tmps = TRUE,
                          seed = 2022,
                          conda_env = "echoR_mini",
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
                          file_sep=deprecated,
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
  # echoverseTemplate:::args2vars(finemap_locus);
  #### Check for required args ####
  force(locus)
  force(fullSS_path)
  ### Check for deprecated args ####
  check_deprecated(fun="finemap_locus",
                   args=match.call())
  #### Create paths ####
  subset_path <- construct_subset_path(results_dir = results_dir,
                                       dataset_type = dataset_type,
                                       dataset_name = dataset_name,
                                       locus = locus)
  locus_dir <- get_locus_dir(subset_path = subset_path)
  ####  Query ####
  steps("query")
  dat <- extract_snp_subset(subset_path = subset_path,
                            locus = locus,
                            topSNPs = topSNPs,
                            fullSS_path = fullSS_path,
                            LD_reference = LD_reference,
                            force_new_subset = force_new_subset,
                            colmap = colmap,
                            bp_distance = bp_distance,
                            superpopulation = superpopulation,
                            compute_n = compute_n,
                            query_by = query_by,
                            nThread = nThread,
                            conda_env = conda_env,
                            verbose = verbose)
  #### Extract LD ####
  steps("ld")
  LD_list <- echoLD::get_LD(query_dat=dat,
                            locus_dir=locus_dir,
                            LD_reference=LD_reference,
                            force_new_LD=force_new_LD,
                            target_genome=LD_genome_build,
                            superpopulation=superpopulation,
                            download_method=download_method,
                            leadSNP_LD_block=leadSNP_LD_block,
                            remove_tmps=remove_tmps,
                            conda_env=conda_env,
                            nThread=nThread,
                            verbose=verbose)
  #### Filter SNPs####
  # Remove pre-specified SNPs
  ## Do this step AFTER saving the LD to disk so that it's easier to re-subset
  ## in different ways later without having to redownload LD.
  LD_list <- echoLD::filter_LD(LD_list=LD_list,
                               remove_correlates=remove_correlates,
                               min_r2=min_r2,
                               verbose=verbose)
  LD_matrix <- LD_list$LD
  dat <- LD_list$DT
  #### Filter SNPs ####
  steps("filter")
  dat <- echodata::filter_snps(
    dat=dat,
    min_MAF=min_MAF,
    bp_distance=bp_distance,
    remove_variants=remove_variants,
    min_POS=min_POS,
    max_POS=max_POS,
    max_snps=max_snps,
    trim_gene_limits=trim_gene_limits,
    verbose=verbose)
  # Subset LD and df to only overlapping SNPs
  LD_list <- echoLD::subset_common_snps(LD_matrix = LD_matrix,
                                        dat = dat,
                                        fillNA = fillNA,
                                        verbose = verbose)
  LD_matrix <- LD_list$LD
  dat <- LD_list$DT
  #### Fine-map ####
  steps("finemap")
  finemap_dat <- echofinemap::multifinemap(locus_dir = locus_dir,
                                 fullSS_path = fullSS_path,
                                 LD_reference = LD_reference,
                                 LD_matrix = LD_matrix,
                                 dat = dat,
                                 # General args (with defaults)
                                 finemap_methods = finemap_methods,
                                 finemap_args = finemap_args,
                                 force_new_finemap = force_new_finemap,
                                 dataset_type = dataset_type,
                                 credset_thresh = credset_thresh,
                                 consensus_thresh = consensus_thresh,
                                 case_control = case_control,
                                 n_causal = n_causal,
                                 compute_n = compute_n,
                                 priors_col = priors_col,
                                 # Tool-specific args
                                 conditioned_snps = conditioned_snps,
                                 # Optional args
                                 seed = seed,
                                 nThread = nThread,
                                 conda_env = conda_env,
                                 verbose = verbose)
  #### LD plot ####
  steps("plot")
  ld_plot <- if("ld" %in% tolower(plot_types)){
    try({
      echoLD::plot_LD(LD_matrix=LD_matrix,
                      dat=finemap_dat,
                      span=10)
    })
  } else {NULL};
  #### Locus plots ####
  locus_plots <- list()
  if("simple" %in% tolower(plot_types)){
    locus_plots[["simple"]] <- tryCatch({
       echoplot::plot_locus(
        dat = finemap_dat,
        LD_matrix = LD_matrix,
        LD_reference = LD_reference,
        locus_dir = locus_dir,
        dataset_type = dataset_type,
        finemap_methods = finemap_methods,
        credset_thresh = credset_thresh,
        consensus_thresh = consensus_thresh,
        tx_biotypes = tx_biotypes,
        qtl_suffixes = qtl_suffixes,
        nott_epigenome = FALSE,
        plot_full_window = TRUE,
        mean.PP = TRUE,
        xgr_libnames = NULL,
        zoom = zoom,
        save_plot = TRUE,
        show_plot = show_plot,
        conda_env = conda_env,
        nThread = nThread,
        verbose = verbose)
    }, error = function(e){message(e);NULL})
  };

  if("fancy" %in% tolower(plot_types)){
    locus_plots[["fancy"]] <- tryCatch({
      echoplot::plot_locus(
        dat = finemap_dat,
        LD_matrix = LD_matrix,
        LD_reference = LD_reference,
        locus_dir = locus_dir,
        dataset_type = dataset_type,
        finemap_methods = finemap_methods,
        credset_thresh = credset_thresh,
        consensus_thresh = consensus_thresh,
        qtl_suffixes = qtl_suffixes,
        zoom = zoom,
        tx_biotypes = tx_biotypes,
        save_plot = TRUE,
        show_plot = show_plot,
        xgr_libnames = xgr_libnames,
        roadmap = roadmap,
        roadmap_query = roadmap_query,
        nott_epigenome = nott_epigenome,
        nott_show_placseq = nott_show_placseq,
        nott_binwidth = nott_binwidth,
        nott_bigwig_dir = nott_bigwig_dir,
        conda_env = conda_env,
        nThread = nThread,
        verbose = verbose)
    }, error = function(e){message(e);NULL})
  };
  #### Return ####
  #### Record args ####
  messager("Recording all `finemap_locus` arguments.",v=verbose)
  # arguments <- tryCatch({
  #   lapply(as.list(match.call(expand.dots=FALSE)),eval)
  # }, error = function(e){message(e);NULL})
  arguments <- NULL
  #### Make list ####
  return(list(finemap_dat=finemap_dat,
              locus_plot=locus_plots,
              LD_matrix=LD_matrix,
              LD_plot=ld_plot,
              locus_dir=locus_dir,
              arguments=arguments
  ))
}

#### Deprecation function #####
finemap_pipeline <- function(...){
  .Deprecated("finemap_locus")
  finemap_locus(...)
}
