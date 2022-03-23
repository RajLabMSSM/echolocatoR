#' Extract a subset of the summary stats
#'
#' Use either \emph{tabix} or \emph{awk} to extract a locus subset
#'  from the full summary statistics file.
#'
#' @family query functions
#' @keywords internal
#' @importFrom echodata check_if_empty
extract_snp_subset <- function(locus=NULL,
                               locus_dir,
                               results_dir=NULL,
                               fullSS_path,
                               fullSS_genome_build="hg19",
                               subset_path,
                               LD_reference,
                               force_new_subset=FALSE,
                               top_SNPs="auto",
                               bp_distance=500000,
                               chrom_col="CHR",
                               chrom_type=NULL,
                               position_col="POS",
                               snp_col="SNP",
                               locus_col="Locus",
                               pval_col="P",
                               effect_col="Effect",
                               stderr_col="StdErr",
                               MAF_col="MAF",
                               freq_col = "Freq",
                               tstat_col="t-stat",
                               A1_col = "A1",
                               A2_col = "A2",
                               gene_col="gene",

                               N_cases_col="N_cases",
                               N_controls_col="N_controls",
                               N_cases=NULL,
                               N_controls=NULL,
                               proportion_cases="calculate",
                               sample_size=NULL,

                               superpopulation="",
                               min_POS=NA,
                               max_POS=NA,
                               genes_detected=FALSE,

                               file_sep="\t",
                               query_by="coordinates",
                               qtl_prefixes=NULL,
                               remove_tmps=TRUE,
                               conda_env = "echoR",
                               verbose=TRUE){
  if(is.null(locus)) locus <- basename(locus_dir)
  multi_path <- echofinemap::create_method_path(
    locus_dir = locus_dir,
    LD_reference = LD_reference,
    finemap_method = "Multi-finemap",
    compress = TRUE)

  if(file.exists(subset_path) & force_new_subset==FALSE){
    messager("+ Importing pre-existing file:",subset_path, v=verbose)
    echodata::check_if_empty(subset_path)
    query <- data.table::fread(subset_path, header=TRUE, stringsAsFactors=FALSE)
  } else if (file.exists(multi_path) & force_new_subset==FALSE){
    messager("+ Importing  pre-existing file:",multi_path, v=verbose)
    echodata::check_if_empty(multi_path)
    query <- data.table::fread(multi_path, header=TRUE, stringsAsFactors=FALSE)
  } else {
    # Extract subset with awk
    messager("+ Extracting relevant variants from fullSS...", v=verbose)
    start_query <- Sys.time()
    # Function selects different methods of querying your SNPs
    query_handler(locus_dir=locus_dir,
                  fullSS_path=fullSS_path,
                  subset_path=subset_path,
                  top_SNPs=top_SNPs,
                  locus_col=locus_col,
                  chrom_col=chrom_col,
                  chrom_type=chrom_type,
                  position_col=position_col,
                  file_sep=file_sep,
                  min_POS=min_POS,
                  max_POS=max_POS,
                  bp_distance=bp_distance,
                  query_by=query_by,
                  force_new_subset=force_new_subset,
                  conda_env = conda_env,
                  verbose=verbose)
    # Clean file
    query <- echodata::standardize_subset(locus=locus,
                                top_SNPs=top_SNPs,
                                fullSS_genome_build=fullSS_genome_build,
                                subset_path=subset_path,
                                chrom_col=chrom_col,
                                position_col=position_col,
                                snp_col=snp_col,
                                pval_col=pval_col,
                                effect_col=effect_col,
                                stderr_col=stderr_col,
                                tstat_col=tstat_col,
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
                                qtl_prefixes=qtl_prefixes)

    end_query <- Sys.time()
    messager("+ Extraction completed in", round(end_query-start_query, 2),
            "seconds", v=verbose)
    messager("+", dim(query)[1], "SNPs x ",dim(query)[2],"columns", v=verbose)
  }
  if(remove_tmps){
    messager("+ Removing subset tmp...", v=verbose)
    suppressWarnings(file.remove(subset_path))
  }
  return(query)
}
