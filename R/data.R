#' Example subset of full summary stats
#'
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @examples
#' A subset of the full GWAS summary stats from Nalls et al. (2019)
#' \dontrun{
#' data("top_SNPs")
#' loci <- c("BST1","LRRK2","MEX3C")
#' fullSS_dat <- lapply(loci, function(locus){
#'   top_sub <- subset(top_SNPs, Locus==locus)
#'   TABIX(fullSS_path = "~/Desktop/nallsEtAl2019_allSamples_allVariants.mod.txt.gz",
#'         save_subset = F,
#'         subset_path = NULL,
#'         chrom_col = "CHR", position_col = "POS",
#'         chrom = top_sub$CHR[1],
#'         min_POS = top_sub$POS - 500000*4,
#'         max_POS = top_sub$POS + 500000*4)
#' }) %>% data.table::rbindlist()
#' usethis::use_data(fullSS_dat, overwrite=T)
#' }
"fullSS_dat"


#' Example subset of full summary stats: munged
#'
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' Munged using \pkg{MungeSumstats}.
#'
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @examples
#' A subset of the full GWAS summary stats from Nalls et al. (2019)
#' \dontrun{
#' data("fullSS_dat")
#' tmp_file <- tempfile(fileext = ".tsv.gz")
#' data.table::fwrite(fullSS_dat, tmp_file)
#' fullSS_munged <- MungeSumstats::format_sumstats(tmp_file, ref_genome = "GRCh37", return_data = TRUE)
#' usethis::use_data(fullSS_munged, overwrite=TRUE)
#' }
"fullSS_munged"



# -------------IMPACT ------------ #


#' IMPACT annotation key
#'
#' Metadata for each of the IMPACT (Inference and Modeling of Phenotype-related ACtive Transcription) annotation files.
#' Originally from \href{https://www.cell.com/ajhg/supplemental/S0002-9297%2819%2930108-9}{Amariuta et al. (2019)}.
#'
#' @family IMPACT
#' @source \url{https://github.com/immunogenomics/IMPACT}
#' @examples
#' \dontrun{
#' IMPACT_annotation_key <- data.table::fread("~/Desktop/Fine_Mapping/echolocatoR/annotations/IMPACT/IMPACT_annotation_key.txt.gz")
#' usethis::use_data(IMPACT_annotation_key, overwrite = T)
#' }
"IMPACT_annotation_key"


