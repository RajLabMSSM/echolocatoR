
#' Column name mappings from \pkg{MungeSumstats} to \pkg{echolocatoR}
#'
#' @source https://github.com/neurogenomics/MungeSumstats
MUNGESUMSTATS.col_map <- function(){
  list(chrom_col="CHR",
    position_col="BP",
    snp_col="SNP",
    pval_col="P",
    effect_col="BETA",
    stderr_col="SE",
    MAF_col="MAF",
    freq_col="FRQ",
    N_cases_col="N_CAS",
    N_controls_col="N_CON",
    A1_col="A1",
    A2_col="A2",
    gene_col="GENE")
}



#' Column name mappings for a given package
#'
#' Supports \pkg{MungeSumstats} and \pkg{echolocatoR}
#'
#' @source https://github.com/neurogenomics/MungeSumstats
column_map <- function(package="MungeSumstats"){
  if(tolower(package)=="mungesumstats"){
    return( list(chrom_col="CHR",
                 position_col="BP",
                 snp_col="SNP",
                 pval_col="P",
                 effect_col="BETA",
                 stderr_col="SE",
                 tstat_col="t_stat",
                 MAF_col="MAF",
                 freq_col="FRQ",
                 N_cases_col="N_CAS",
                 N_controls_col="N_CON",
                 A1_col="A1",
                 A2_col="A2"))
  }
  if(tolower(package)=="echolocator"){
    return( list(chrom_col="CHR",
                 position_col="POS",
                 snp_col="SNP",
                 pval_col="P",
                 effect_col="Effect",
                 stderr_col="StdErr",
                 tstat_col="t_stat",
                 MAF_col="MAF",
                 freq_col="Freq",
                 N_cases_col="N_cases",
                 N_controls_col="N_controls",
                 A1_col="A1",
                 A2_col="A2",
                 locus_col="Locus",
                 gene_col="Gene"))
  }

}


#' Convert from \pkg{MungeSumstats} to \pkg{echolocatoR} format
#'
#' @examples
#' \dontrun{
#' eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt", package="MungeSumstats")
#' reformatted <- MungeSumstats::format_sumstats(path=eduAttainOkbayPth, ref_genome="GRCh37", return_data = TRUE)
#' dat_echoR <- MUNGESUMSTATS.to_echolocatoR(dat=reformatted)
#' }
MUNGESUMSTATS.to_echolocatoR <- function(dat){
  printer("+ Mapping colnames from MungeSumstats ==> echolocatoR")
  dat <- data.table::data.table(dat)
  colMap <- c("BP"="POS","FRQ"="Freq","BETA"="Effect","SE"="StdErr","N_CAS"="N_cases","N_CON"="N_controls")
  colMap <- colMap[names(colMap) %in% colnames(dat)]
  for(x in names(colMap)){
    setnames(dat, x, colMap[[x]])
  }
  return(dat)
}



#' Check whether a table has headers that can be mapped to an attribute in \pkg{MungeSumstats}
#'
#' @family internal
MUNGESUMSTATS.check_syn <- function(dat,
                                    col_name="P-value",
                                    col_type="P"){
  if(any(col_name==col_type)) return(dat)
  if(!is.null(col_name)){
     dat <- data.frame(dat, check.names = FALSE)
    if(col_name!=col_type){
      dat[[col_type]] <- dat[[col_name]]
    }
  }
  syn <- subset(MungeSumstats:::sumstatsColHeaders, Corrected==col_type)[,1]
  if(!any(colnames(dat) %in% syn)) stop(paste("MungeSumstats will be unable to map",col_type,", please supply mapping."))
  return(dat)
}
