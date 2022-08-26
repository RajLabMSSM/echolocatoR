#' Convert Jack's coloc results to \emph{echolocatoR} format
#'
#' @param all_obj Nested list object created by Jack.
#' @param results_level Return coloc results at the
#'  "summary" (one row per Locus:eGene pair) or "snp" level (one row per SNP).
#' @param save_path File where you want the merged results to be saved.
#' @family coloc
#' @keywords internal
#' @examples
#' \dontrun{
#' load("/sc/hydra/projects/ad-omics/microglia_omics/COLOC/Kunkle_2019/Microglia_all_regions_Kunkle_2019_COLOC.RData")
#' merged_results <- merge_coloc_results(all_obj=all_obj, results_level="snp", save_path="~/Desktop")
#' }
merge_coloc_results <- function(all_obj,
                                results_level=c("summary"),
                                nThread=1,
                                verbose=TRUE,
                                save_path=FALSE){
  null_list<<-NULL
  messager("Gathering coloc results at",results_level,"level...")
  merged_results <- parallel::mclapply(names(all_obj), function(locus){
    messager("- Locus =",locus,v=verbose)
    locus_obj <- all_obj[[locus]]
    parallel::mclapply(names(locus_obj), function(egene){
      messager("eGene =",egene,v=verbose)
      egene_obj <- locus_obj[[egene]]
      if(results_level=="snp"){
        results <- egene_obj$object$results
        if(!is.null(results)){
          results <- cbind(Locus=locus,results)
        } else {null_list <<- append(null_list, setNames(locus, egene))}
      } else{
        results <- egene_obj$df
        if(!is.null(results)){
          results <- cbind(Locus=locus,gene=egene,results)
        } else {null_list <<- append(null_list, setNames(locus, egene))}
      }
      return(results)
    }, mc.cores=nThread) %>% data.table::rbindlist()
  }, mc.cores = 1)  %>% data.table::rbindlist()
  if(length(null_list)>0) {
    messager("NULL results detected in",length(null_list),"Locus:eGene pairs.");
    print(null_list);
  }
  # Rename cols
  suffix_to_prefix <- function(dat, suffix){
    cols_select <- grep(paste0(suffix,"$"),colnames(dat))
    prefix <- paste0(gsub("\\.","",suffix),".")
    colnames(dat)[cols_select] <- paste0(prefix, gsub(paste0(suffix,"$"),"",colnames(dat)[cols_select] ))
    return(dat)
  }
  merged_results <- suffix_to_prefix(dat=merged_results, suffix = ".gwas")
  merged_results <- suffix_to_prefix(dat=merged_results, suffix = ".qtl")

  # get N_cases and N_controls from N.gwas and proportion_controls
  merged_results$N_cases <- floor(merged_results$gwas.N * merged_results$s)
  merged_results$N_controls <- merged_results$gwas.N - merged_results$N_cases

  if(save_path!=FALSE){
    data.table::fwrite(merged_results,
                       save_path,
                       sep="\t", nThread=nThread)
  }
  return(merged_results)
}

