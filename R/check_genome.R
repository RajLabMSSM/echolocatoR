check_genome <- function(gbuild=NULL,
                         munged=FALSE,
                         fullSS_path=NULL){
  if(is.null(gbuild)){
    if(munged & (!is.null(fullSS_path))){
      messager("+ Inferring genome build")
      sumstats <- data.table::fread(fullSS_path, nThread = 1, nrows = 1000)
      gbuild <- MungeSumstats:::get_genome_build(sumstats = sumstats)
      return(gbuild)
    } else {
      message("WARNING:: fullSS_genome_build not provided. Assuming hg19.")
      gbuild <- "GRCH37"
      return(gbuild)
    }
  } else {
    gbuild <- if(tolower(gbuild) %in% tolower(c("hg19","hg37","GRCh37","grch37"))) "GRCH37" else gbuild
    gbuild <- if(tolower(gbuild) %in% tolower(c("hg38","GRCh38","grch38"))) "GRCH38" else gbuild
    return(gbuild)
  }
}
