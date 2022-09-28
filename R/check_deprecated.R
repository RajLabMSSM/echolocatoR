#' Check deprecated arguments
#'
#' Semi-automatically check all deprecated args in a given function.
#' @param fun Function to check.
#' @param pkg Package that the function is from.
#' @param args Argument calls to assess.
#' @param lifecycle_fun Which \pkg{lifecycle} function to use by default.
#' @param reassign Attempt to reassign deprecated variables to
#' the corresponding new variable (if applicable).
#' @param map Mapping between old:new argument names. Use \code{NULL}
#' if the argument is no longer used at all.
#' @inheritParams lifecycle::deprecate_warn
#'
#' @export
#' @importFrom stats na.omit
#' @importFrom lifecycle deprecate_stop
#' @importFrom utils getFromNamespace
#' @examples
#' topSNPs <- echodata::topSNPs_Nalls2019
#' fullSS_path <- echodata::example_fullSS()
#' testthat::expect_error(
#'   echolocatoR::finemap_loci(
#'     fullSS_path = fullSS_path,
#'     topSNPs = topSNPs,
#'     loci = c("BST1","MEX3C"),
#'     chrom_col = "CHR",
#'     position_col = "BP")
#' )
check_deprecated <- function(fun="finemap_loci",
                             pkg="echolocatoR",
                             when="2.0.0",
                             args=match.call(),
                             lifecycle_fun=lifecycle::deprecate_warn,
                             reassign=FALSE,
                             map=list(A1_col="colmap",
                                     A2_col="colmap",
                                     chrom_col="colmap",
                                     position_col="colmap",
                                     effect_col="colmap",
                                     freq_col="colmap",
                                     gene_col="colmap",
                                     locus_col="colmap",
                                     MAF_col="colmap",
                                     N_cases="colmap",
                                     N_controls="colmap",
                                     N_cases_col="colmap",
                                     N_controls_col="colmap",
                                     sample_size="colmap",
                                     MAF_col="colmap",
                                     pval_col="colmap",
                                     stderr_col="colmap",
                                     tstat_col="colmap",
                                     snp_col="colmap",
                                     file_sep=NULL,
                                     probe_path=NULL,
                                     chrom_type=NULL,
                                     PAINTOR_QTL_datasets=NULL,
                                     QTL_prefixes="qtl_suffixes",
                                     proportion_cases=NULL,
                                     server=NULL,
                                     vcf_folder=NULL,
                                     top_SNPs="topSNPs",
                                     PP_threshold="credset_thresh",
                                     consensus_threshold="consensus_thresh",
                                     plot.types="plot_types",
                                     plot.Roadmap="roadmap",
                                     plot.Roadmap_query="roadmap_query",
                                     plot.XGR_libnames="xgr_libnames",
                                     plot.zoom="zoom",
                                     plot.zoom="zoom",
                                     plot.Nott_epigenome="nott_epigenome",
                                     plot.Nott_show_placseq="nott_show_placseq"
                                     )
                             ){

  #### Test run ####
  # echoverseTemplate:::args2vars(check_deprecated)
  # fun <- function(x, lower = 0, upper = 1) {
  #   as.list(match.call())
  # }
  # args <- fun(4 * atan(1), u = pi)
  # names(args)[-1] <- c("PAINTOR_QTL_datasets","A2_col")#names(map)[seq_len(2)]


  #### Find all deprecated args ####
  dep_args <- formals(utils::getFromNamespace(fun,pkg))
  dep_args <- sort(
    as.character(
      stats::na.omit(
        names(
          unlist(dep_args)[
            unlist(lapply(dep_args, as.character))=="deprecated"
            ]
        )
      )
    )
  )
  args <- as.list(args)[-1]
  args <- args[names(args) %in% dep_args]
  if(length(args)>0){
    wng <- paste("Found",length(args),"deprecated arguments used in:",
                 deparse(args[[1]]))
    warning(wng)
    for(nm in names(args)){
      if (lifecycle::is_present(args[[nm]])) {
        ##### Check with ####
        with <- if(nm %in% names(map)) {
          paste0(pkg,"::",fun,"(",map[[nm]],")")
        } else {
          paste0(pkg,"::",fun,"()")
        }
        #### Run appropriate warning/error ####
        if(is.null(map[[nm]])){
          lifecycle_fun(when = when,
                        what = paste0(fun,"(",nm,")"))
        }else if(map[[nm]]=="colmap"){
            lifecycle::deprecate_stop(when = when,
                                      what = paste0(fun,"(",nm,")"),
                                      with = with)

        }else{
          lifecycle_fun(when = when,
                        what = paste0(fun,"(",nm,")"),
                        with = with)
          #### Reassign values to variables ####
          if(isTRUE(reassign)){
            messager("Reassigning deprecated argument:",nm,"-->",map[[nm]])
            assign(x = map[[nm]],
                   value = get(nm),
                   inherits = TRUE,
                   pos = 1L)
          } else {
            messager("Deprecated argument will be ignored:",nm)
          }
        }
      }
    }
  }
}
