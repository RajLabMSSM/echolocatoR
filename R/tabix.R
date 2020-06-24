


#' Convert summary stats file to tabix format
#'
#' @family query functions
#' @inheritParams finemap_pipeline
#' @examples
#' \dontrun{
#' data("genome_wide_dir");
#' fullSS_path <-  "~/Desktop/nallsEtAl2019_allSamples_allVariants.mod.txt"
#' fullSS.gz <- TABIX.convert_file(fullSS_path=fullSS_path, results_path=results_path_genome_wide, chrom_col="CHR", position_col="POS")
#' }
TABIX.convert_file <- function(fullSS_path,
                               results_path=NULL,
                               chrom_col="CHR",
                               position_col="POS",
                               verbose=T){
  results_path <- ifelse(is.null(results_path),dirname(fullSS_path),results_path)
  fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), file.path(results_path, basename(fullSS_path)), paste0(file.path(results_path, basename(fullSS_path)),".gz"))
  z_grep <- ifelse(endsWith(fullSS_path,".gz"),"zgrep","grep")
  cDict <-  column_dictionary(file_path = fullSS_path) # header.path
  print("TABIX:: Converting full summary stats file to tabix format for fast querying...")
  cmd <- paste("(",
               # Extract the header col and sort everything else
               paste0(z_grep," '",chrom_col,"' ",fullSS_path,"; ",z_grep," -v ^'",chrom_col,"' ",fullSS_path),
               paste0("| sort -k",cDict[[chrom_col]],",",cDict[[chrom_col]]),
               paste0("-k",cDict[[position_col]],",",cDict[[position_col]],"n"),
               ")",
               # Compress with bgzip
               "| bgzip >",
               fullSS.gz)
  printer(cmd, v=verbose)
  system(cmd)

  # Index
  TABIX.index <- function(fullSS.gz,
                          chrom_i=1,
                          pos_i=2,
                          skip_lines=1,
                          verbose=T){
    printer("TABIX:: Indexing ")
    cmd2 <- paste("tabix",
                  "-f", # Force overwrite of .tbi index file
                  "-S",skip_lines,#--skip-lines
                  "-s",chrom_i,
                  "-b",pos_i,
                  "-e",pos_i,
                  fullSS.gz)
    printer(cmd2, v=verbose)
    system(cmd2)
  }
  TABIX.index(fullSS.gz=fullSS.gz,
              chrom_i=cDict[[chrom_col]],
              pos_i=cDict[[position_col]],
              skip_lines = 1,
              verbose=verbose)
  return(fullSS.gz)
}


#' Query a tabix file
#'
#' Query by genomic coordinates.
#'
#' @family query functions
TABIX.query <- function(fullSS.gz,
                        chrom,
                        start_pos,
                        end_pos,
                        verbose=T){
  coords <- paste0(chrom,":",start_pos,"-",end_pos)
  # cmd4 <- paste("tabix -h",fullSS.gz,coords,">",subset_path)
  printer("TABIX:: Extracting subset of sum stats", v=verbose)
  dat <- data.table::fread(cmd=paste("tabix -h",fullSS.gz,coords))
  printer("++ Returning",paste(dim(dat),collapse=" x "),"data.table", v=verbose)
  return(dat)
}




#' Covert and query
#'
#' If it is not tabix format already
#' (determined by checking for a \link{.tbi} file of the same name in the same directory),
#' the full summary statistics file is converted into tabix format for super fast querying.
#' A query is then made using the min/max genomic positions to extract a locus-specific summary stats file.
#'
#' @family query functions
#' @inheritParams finemap_pipeline
#' @return data.table of locus subset summary statistics
#' @examples
#' \dontrun{
#' data("locus_dir"); data("Nalls_top_SNPs")
#' fullSS_path <- "./results/GWAS/Nalls23andMe_2019/_genome_wide/nallsEtAl2019_allSamples_allVariants.mod.txt.gz"
#' top_SNPs <- import_topSNPs(topSS = Nalls_top_SNPs, chrom_col = "CHR", position_col = "BP", snp_col="SNP", pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene", group_by_locus = T,locus_col = "Nearest Gene")
#' top_SNPs_BST1 <- subset(top_SNPs, Locus=='BST1')
#' bp_distance <- 1e+06
#' min_POS <- top_SNPs_BST1$POS - bp_distance
#' max_POS <- top_SNPs_BST1$POS + bp_distance
#' subset_file <- file.path(results_path,"BST1_Nalls23andMe_2019_subset.tsv.gz")
#' dat <- TABIX(fullSS_path=fullSS_path, subset_path="auto", min_POS=min_POS, max_POS=max_POS, chrom=top_SNPs_BST1$CHR)
#' }
TABIX <- function(fullSS_path,
                  subset_path=NULL,
                  is_tabix=F,
                  chrom_col="CHR",
                  position_col="POS",
                  min_POS=NA,
                  max_POS=NA,
                  chrom=NULL,
                  save_subset=T,
                  nThread=1,
                  verbose=T){
  # Check if it's already an indexed tabix file
  fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), fullSS_path, paste0(fullSS_path,".gz"))
  is_tabix <- file.exists(fullSS.gz) & file.exists(paste0(fullSS.gz,".tbi"))
  if(!is_tabix){
    fullSS.gz <- TABIX.convert_file(fullSS_path,
                                    chrom_col = chrom_col,
                                    position_col = position_col,
                                    verbose=verbose)
  } else { printer("TABIX:: Existing indexed tabix file detected") }
  # Query
  cDict <- column_dictionary(file_path = fullSS.gz)
  dat <- TABIX.query(fullSS.gz=fullSS.gz,
                      chrom=chrom,
                      start_pos=min_POS,
                      end_pos=max_POS,
                      verbose=verbose)
  colnames(dat) <- names(cDict)
  if(save_subset){
    printer("++ Saving query ==>", subset_path, v=verbose)
    data.table::fwrite(dat, file = subset_path, nThread = nThread, sep="\t")
  }
  return(dat)
}




