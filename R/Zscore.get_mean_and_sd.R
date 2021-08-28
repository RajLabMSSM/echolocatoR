#' Store info necessary for Z-score
#'
#' Store the info about the full vector (length, mean, standard deviation),
#' so that you can later accurately compute Z-scores
#' when you're only working with a subset of the data.
#'
#' These functions are necessary for \code{\link{PAINTOR}}.
#' @keywords internal
Zscore.get_mean_and_sd <- function(fullSS,#="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                                   target_col="statistic",
                                   effect_col="beta",
                                   stderr_col="se",
                                   use_saved=T,
                                   output_path){#="./Data/GWAS/Nalls23andMe_2019/z.info.RDS"){
  if(use_saved & file.exists(output_path)){
    printer("Reading in:",output_path,"...")
    z.info <- readRDS(output_path)
  } else {
    printer("Extracting mean and standard deviation from",fullSS,"...")
    if(target_col=="calculate"){
      target_col <- "t_stat"
      sample_x <- data.table::fread(fullSS, nThread = 4,
                                    select=c(effect_col, stderr_col),
                                    col.names = c("Effect","StdErr"))
      sample_x <- subset(calculate_tstat(sample_x), select = target_col)
    }else {
      sample_x <- data.table::fread(fullSS, nThread = 4, select=c(target_col))
    }
    sample.mean <- mean(sample_x[[1]], na.rm = T)
    sample.stdv <- sd(sample_x[[1]])
    z.info <- list(file.name=fullSS,
                   colname=target_col,
                   sample.mean=sample.mean,
                   sample.stdv=sample.stdv,
                   sample.min=min(sample_x[[1]]),
                   sample.max=max(sample_x[[1]]))
    saveRDS(z.info, file = output_path )
  }
  return(z.info)
}
