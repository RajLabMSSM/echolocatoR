#' IMPACT annotation key
#'
#' Metadata for each of the IMPACT
#' (Inference and Modeling of Phenotype-related ACtive Transcription)
#' annotation files.
#' Originally from
#' \href{https://www.cell.com/ajhg/supplemental/S0002-9297%2819%2930108-9}{
#' Amariuta et al. (2019)}.
#'
#' @family IMPACT
#' @source \url{https://github.com/immunogenomics/IMPACT}
#' @source
#' \code{
#' path <- file.path("~/Desktop/Fine_Mapping/echolocatoR",
#'                   "annotations/IMPACT/IMPACT_annotation_key.txt.gz")
#' IMPACT_annotation_key <- data.table::fread(path)
#' usethis::use_data(IMPACT_annotation_key, overwrite = TRUE)
#' }
"IMPACT_annotation_key"
