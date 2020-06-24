#' @details
#' Fine-mapping methods are a powerful means of identifying causal variants underlying a given phenotype,
#' but are  underutilized due to the technical challenges of implementation.
#' \emph{echolocatoR} is an R package that automates end-to-end genomics fine-mapping, annotation,
#' and plotting in order to identify the most probable causal variants associated with a given phenotype.
#'
#' It requires minimal input from users (a GWAS or QTL summary statistics file),
#' and includes a suite of statistical and functional fine-mapping tools.
#' It also includes extensive access to datasets
#' (linkage disequilibrium panels, epigenomic and genome-wide annotations, QTL).
#'
#' The elimination of data gathering and preprocessing steps enables rapid fine-mapping of many loci in any phenotype,
#'  complete with locus-specific publication-ready figure generation.
#'  All results are merged into a single per-SNP summary file for additional downstream analysis
#'   and results sharing. Therefore \emph{echolocatoR} drastically reduces the barriers to identifying
#'   causal variants by making the entire fine-mapping pipeline rapid, robust and scalable.
#' @keywords internal
"_PACKAGE"
