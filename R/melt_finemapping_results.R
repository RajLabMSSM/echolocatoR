melt_finemapping_results <- function(finemap_dat,
                                     verbose=T){
  PP_cols <- grep("\\.PP$",colnames(finemap_dat), value = T)
  CS_cols <- grep("\\.CS$",colnames(finemap_dat), value = T)
  printer("Melting PP and CS from",length(CS_cols),"fine-mapping methods",v=verbose)
  finemap_melt <- suppressWarnings(
    data.table::melt.data.table(data.table::data.table(finemap_dat),
                                measure.vars =list(PP_cols, CS_cols),
                                variable.name = "Method",
                                variable.factor = T,
                                value.name = c("PP","CS"))
  )
  methods_key <- setNames(gsub("\\.PP","",PP_cols),
                          unique(finemap_melt$Method))
  finemap_melt$Method <- factor(methods_key[finemap_melt$Method], levels = unname(methods_key), ordered = T)
  return(finemap_melt)
}
