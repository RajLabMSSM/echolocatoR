#' Run \pkg{echolocatoR} pipeline on a single locus
#'
#' Unlike \code{finemap_loci}, you don't need to provide a \code{top_SNPs}
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
#' @section input file parameters:
#'
#' @param loci Character list of loci in \strong{Locus} col of \code{top_SNPs}.
#' @param fullSS_path Path to the full summary statistics file (GWAS or QTL)
#' that you want to fine-map.
#' It is usually best to provide the absolute path rather
#' than the relative path.
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
#' @param top_SNPs A data.frame with the genomic coordinates of the lead SNP
#' for each locus.
#' The lead SNP will be used as the center of the window when extracting
#' subset from the full GWAS/QTL summary statistics file.
#' Only one SNP per \strong{Locus} should be included.
#' At minimum, \code{top_SNPs} should include the following columns:
#' \describe{
#' \item{\emph{Locus}}{A unique name for each locus. Often,
#'  loci are named after a relevant gene (e.g. LRRK2) or based on
#'   the name/coordinates of the lead SNP (e.g. locus_chr12_40734202) }
#' \item{\emph{CHR}}{The chromosome that the SNP is on.
#'  Can be "chr12" or "12" format.}
#' \item{\emph{POS}}{The genomic position of the SNP (in basepairs)}
#' }
#'
#' @section overwrite existing files:
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
#' @section fine-mapping parameters:
#'
#' @param finemap_methods Which fine-mapping methods you want to use.
#' @param n_causal The maximum number of potential causal SNPs per locus.
#' This parameter is used somewhat differntly by different fine-mapping tools.
#' See tool-specific functions for details.
#' @param munged Whether \code{fullSS_path} have already been
#' standardised/filtered  full summary stats
#' with \link[MungeSumstats]{format_sumstats}.
#' If \code{munged=FALSE} you'll need to provide the necessary
#'  column name arguments.
#' @param conditioned_snps Which SNPs to conditions on when fine-mapping
#' with \emph{COJO}.
#' @param PAINTOR_QTL_datasets A list of QTL datasets to be used when
#' conducting joint functional fine-mapping with \emph{PAINTOR}.
#' @param PP_threshold The minimum fine-mapped posterior probability
#'  for a SNP to be considered part of a Credible Set.
#' For example, \code{PP_threshold=.95} means that all Credible Set SNPs
#' will be 95\% Credible Set SNPs.
#' @param consensus_threshold The minimum number of fine-mapping tools
#' that include a SNP
#'  in their 95\% Credible Sets to consider that it a "Consensus SNP"
#'  (\emph{default=2}).
#'
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
#'
#' @family MAIN
#' @inheritParams echoconda::activate_env
#' @inheritParams echodata::filter_snps
#' @inheritParams echoLD::get_LD
#' @inheritParams echoLD::filter_LD
#' @inheritParams echoplot::plot_locus
#' @inheritParams echofinemap::multifinemap
#' @importFrom echoplot plot_locus
#' @importFrom echoLD get_LD filter_LD subset_common_snps
#' @importFrom echofinemap multifinemap
#' @export
#' @examples
#' top_SNPs <- echodata::topSNPs_Nalls2019
#' fullSS_path <- echodata::example_fullSS(dataset = "Nalls2019")
#'
#' res <- echolocatoR::finemap_locus(
#'   fullSS_path = fullSS_path,
#'   top_SNPs = top_SNPs,
#'   locus = "BST1",
#'   dataset_name = "Nalls23andMe_2019",
#'   fullSS_genome_build = "hg19",
#'   bp_distance=10000,
#'   munged = TRUE)
finemap_locus <- function(#### Main args ####
                          locus,
                          fullSS_path,
                          fullSS_genome_build=NULL,
                          results_dir=file.path(tempdir(),"results"),
                          dataset_name="dataset_name",
                          dataset_type="GWAS",
                          case_control=TRUE,
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
                          conditioned_snps=NULL,
                          #### Colname mapping args ####
                          munged = FALSE,
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
                          show_plot=TRUE,
                          #### General args ####
                          remove_tmps=TRUE,
                          conda_env="echoR",
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
                          sample_size=deprecated()
                          ){
  #### Map columns from MungeSumstats ####
  if(munged){
    messager(
      "+ Assuming sumstats have already been processed with MungeSumstats.",
      v=verbose)
  }
  #### Create paths ####
  subset_path <- get_subset_path(results_dir = results_dir,
                                 dataset_type = dataset_type,
                                 dataset_name = dataset_name,
                                 locus = locus)
  locus_dir <- get_locus_dir(subset_path = subset_path)

  ####  Query ####
  dat <- extract_snp_subset(locus = locus,
                                  results_dir = results_dir,
                                  locus_dir = locus_dir,
                                  top_SNPs = top_SNPs,
                                  fullSS_path = fullSS_path,
                                  fullSS_genome_build = fullSS_genome_build,
                                  subset_path  =  subset_path,
                                  LD_reference = LD_reference,
                                  force_new_subset = force_new_subset,

                                  chrom_col = chrom_col,
                                  chrom_type = chrom_type,
                                  position_col = position_col,
                                  snp_col = snp_col,
                                  pval_col = pval_col,
                                  effect_col = effect_col,
                                  stderr_col = stderr_col,
                                  locus_col = locus_col,
                                  tstat_col = tstat_col,
                                  MAF_col = MAF_col,
                                  freq_col = freq_col,
                                  A1_col = A1_col,
                                  A2_col = A2_col,
                                  gene_col = gene_col,

                                  N_cases_col = N_cases_col,
                                  N_controls_col = N_controls_col,
                                  N_cases = N_cases,
                                  N_controls = N_controls,
                                  proportion_cases = proportion_cases,
                                  sample_size = sample_size,

                                  bp_distance = bp_distance,
                                  superpopulation = superpopulation,
                                  min_POS = min_POS,
                                  max_POS = max_POS,

                                  query_by = query_by,
                                  qtl_prefixes = qtl_prefixes,

                                  remove_tmps = remove_tmps,
                                  conda_env = conda_env,
                                  verbose = verbose)
  #### Extract LD ####
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
  finemap_dat <- echofinemap::multifinemap(locus_dir = locus_dir,
                                 fullSS_path = fullSS_path,
                                 LD_reference = LD_reference,
                                 LD_matrix = LD_matrix,
                                 dat = dat,
                                 sample_size = sample_size,
                                 # General args (with defaults)
                                 finemap_methods = finemap_methods,
                                 finemap_args = finemap_args,
                                 force_new_finemap = force_new_finemap,
                                 dataset_type = dataset_type,
                                 PP_threshold = PP_threshold,
                                 consensus_threshold = consensus_threshold,
                                 case_control = case_control,
                                 n_causal = n_causal,
                                 # Tool-specific args
                                 conditioned_snps = conditioned_snps,
                                 PAINTOR_QTL_datasets = PAINTOR_QTL_datasets,
                                 # Optional args
                                 conda_env = conda_env,
                                 verbose = verbose)
  #### LD plot ####
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
    try({
      TRKS_simple <- echoplot::plot_locus(
        dat = finemap_dat,
        LD_matrix = LD_matrix,
        LD_reference = LD_reference,
        locus_dir = locus_dir,
        dataset_type = dataset_type,
        finemap_methods = finemap_methods,
        PP_threshold = PP_threshold,
        consensus_threshold = consensus_threshold,
        qtl_prefixes = qtl_prefixes,
        nott_epigenome = FALSE,
        plot_full_window = FALSE,
        mean.PP = TRUE,
        xgr_libnames = NULL,
        zoom = zoom,
        save_plot = TRUE,
        show_plot = show_plot,
        conda_env = conda_env,
        nThread = nThread,
        verbose = verbose)
      locus_plots[["simple"]] <- TRKS_simple
    })
  };

  if("fancy" %in% tolower(plot_types)){
    try({
      TRKS_fancy <- echoplot::plot_locus(
        dat = finemap_dat,
        LD_matrix = LD_matrix,
        LD_reference = LD_reference,
        locus_dir = locus_dir,
        dataset_type = dataset_type,
        finemap_methods = finemap_methods,
        PP_threshold = PP_threshold,
        consensus_threshold = consensus_threshold,
        qtl_prefixes = qtl_prefixes,
        zoom = zoom,
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
      locus_plots[["fancy"]] <- TRKS_fancy
    })
  };
  #### Cleanup ####
  if(remove_tmps){
    roadmap_tbi <- list.files(locus_dir,pattern=".bgz.tbi",
                              recursive = TRUE, full.names = TRUE)
    tmp_files <- file.path(locus_dir,"LD",
                           c("plink.bed",
                             "plink.bim",
                             "plink.fam",
                             "plink.ld",
                             "plink.ld.bin",
                             "plink.log",
                             "plink.nosex",
                             "SNPs.txt") )
    out <- suppressWarnings(file.remove(tmp_files, roadmap_tbi))
  }
  return(list(finemap_dat=finemap_dat,
              locus_plot=locus_plots,
              LD_matrix=LD_matrix,
              LD_plot=ld_plot,
              locus_dir=locus_dir,
              arguments=as.list(match.call(expand.dots=FALSE))
  ))
}

#### Deprecation function #####
finemap_pipeline <- function(...){
  .Deprecated("finemap_locus")
  finemap_locus(...)
}
