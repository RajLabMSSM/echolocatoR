#' Impute t-stat
#'
#' If \strong{tstat} column is missing,
#' compute t-statistic from: \code{Effect / StdErr}.
#' @family standardization functions
calculate_tstat <- function(finemap_dat,
                            tstat_col="t_stat"){
  if(tstat_col %in% colnames(finemap_dat)){
    finemap_dat <- finemap_dat %>% dplyr::rename(t_stat = tstat_col)
  } else {
    if(("Effect" %in% colnames(finemap_dat)) &
       ("StdErr" %in% colnames(finemap_dat))){
      messager("+ Calculating t-statistic from Effect and StdErr...")
      finemap_dat <- finemap_dat %>% dplyr::mutate(t_stat =  Effect/StdErr)
    } else {
      messager(
        "+ Could not calculate t-stat",
        "due to missing Effect and/or StdErr columns. Returning input data.")
    }
  }
  return(data.table::data.table(finemap_dat))
}
