

# BiocManager::install("motifbreakR")
# BiocManager::install(c("SNPlocs.Hsapiens.dbSNP142.GRCh37","BSgenome.Hsapiens.UCSC.hg19"))
# library(motifbreakR); library(BSgenome)


#' Run \code{\link{motifbreakR}}
#'
#' \code{\link{motifbreakR}} is a package to predict how much a SNP will disrupt
#' a transcription factor binding motif (if it falls witihn one).
#' @return motifbreakr results
#' @source
#' \strong{Publication:}
#' \url{https://pubmed.ncbi.nlm.nih.gov/26272984/}
#' \strong{GitHub:}
#' \url{https://github.com/Simon-Coetzee/MotifBreakR}
#' \strong{Vignette:}
#' \url{http://simon-coetzee.github.io/motifBreakR}
#' @examples
#' \dontrun{
#' data("merged_DT")
#' library('BSgenome')
#' # merged_DT <- openxlsx::read.xlsx("/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/Nalls23andMe_2019_results.xlsx")
#' snp_list <- unique(merged_DT$SNP)
#' motifbreakr.results <- MOTIFBREAKR(snp_list=snp_list)
#' }
#'
MOTIFBREAKR <- function(snp_list,
                        save_rds=T,
                        dataset_dir){
  # save_rds=T; dataset_dir <- "./results/GWAS/Nalls23andMe_2019"; library(BSgenome);

  variants <- motifbreakR::snps.from.rsid(rsid = snp_list,
                                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                          search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19);
  motifbreakr.results <- motifbreakR::motifbreakR(snpList = variants,
                                                  pwmList = MotifDb::MotifDb,
                                                  threshold = 0.9);
  if(save_rds){
    rds_path <- file.path(dataset_dir,'_genome_wide','motifbreakR','motifbreakR_results.rds');
    dir.create(dirname(rds_path),showWarnings = F, recursive = T);
    # printer("+ MOTIFBREAKR:: Saving results ==>", rds_path);
    saveRDS(motifbreakr.results, rds_path);
  }
  return(motifbreakr.results)
}




#' Plot \code{\link{motifbreakR}} results
#'
#' @source
#' \strong{Publication:}
#' \url{https://pubmed.ncbi.nlm.nih.gov/26272984/}
#'
#' \strong{GitHub:}
#' \url{https://github.com/Simon-Coetzee/MotifBreakR}
#'
#' @examples
#' \dontrun{
#' # motifbreakr.results <- readRDS("/sc/arion/projects/pd-omics/brian/motifbreakR/motifbreakR_results.rds")
#' motifbreakr.results <- readRDS("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/motifbreakR/motifbreakR_results.rds")
#' MOTIFBREAKR.plot(motifbreakr.results=motifbreakr.results, rsid="rs114528427")
#' }
MOTIFBREAKR.plot <- function(motifbreakr.results,
                             rsid=NULL){
  library(BSgenome); library(BSgenome.Hsapiens.UCSC.hg19); library(motifbreakR);
  # rsid<-"rs7294619"
  if(is.null(rsid)){rsid <- names(motifbreakr.results)[5]}

  motifbreakR::plotMB(results = motifbreakr.results,
                      rsid = rsid,
                      effect = "strong")
}

