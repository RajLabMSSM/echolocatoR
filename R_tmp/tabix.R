


TABIX.convert_file <- function(fullSS_path="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                               chrom_col="CHR",
                               position_col="POS",
                               remove_header=T){
  # Instantaneous header collection without reading in large file!
  # headers <- colnames(data.table::fread(cmd=paste("head -1",fullSS_path)))
  header.path <- file.path(dirname(fullSS_path),"header.txt")
  fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), fullSS_path, paste0(fullSS_path,".gz"))
  gz_cat <- ifelse(endsWith(fullSS_path,".gz"),"gunzip -c","cat")
  # Get header and cdict
  if(!file.exists(header.path)){
    system( paste(gz_cat,fullSS_path, "| head -1 >",header.path) )
  }
  cDict <-  column_dictionary(file_path = header.path) # header.path
  print("TABIX:: Converting full summary stats file to tabix format for fast querying...")

  cmd <- paste(gz_cat,
               fullSS_path,
               # Get rid of header
               ifelse(remove_header,"| tail -n +2",""),
               # Sort (sort -k1,1 -k2,2n)
               paste0("| sort -k",cDict[[chrom_col]],",",cDict[[chrom_col]]),
               paste0("-k",cDict[[position_col]],",",cDict[[position_col]],"n"),
               # ">",fullSS_path)
               # Compress with bgzip
               "| bgzip >",
               # Must save a tmp file first, and them rename/move to final destination
               " tmp && mv tmp",fullSS.gz)
  print(cmd)
  system(cmd)
  # Rsamtools:::bgzip(file = fullSS_path,
  #                   dest = fullSS.gz,
  #                   overwrite = T)
  # idx <- Rsamtools::indexTabix(file = fullSS.gz,
  #                              format = "bed",
  #                              seq = 1,
  #                              start = 2,
  #                              end = 2,
  #                              skip=1)
  # Rsamtools::headerTabix(tbx)
  # tab <- TabixFile(zipped, idx)


  # Index
  TABIX.index <- function(fullSS.gz,
                          chrom_i=1,
                          pos_i=2,
                          skip_lines=1){
    printer("TABIX:: Indexing ")
    cmd2 <- paste("tabix",
                  "-f", # Force overwrite of .tbi index file
                  "-S",skip_lines,#--skip-lines
                  "-s",chrom_i,
                  "-b",pos_i,
                  "-e",pos_i,
                  fullSS.gz)
    print(cmd2)
    system(cmd2)
  }
  TABIX.index(fullSS.gz=fullSS.gz,
              chrom_i=cDict[[chrom_col]],
              pos_i=cDict[[position_col]],
              skip_lines = 1)
  return(fullSS.gz)
}


# Query
TABIX.query <- function(fullSS.gz,
                        chrom,
                        start_pos,
                        end_pos){
  coords <- paste0(chrom,":",start_pos,"-",end_pos)
  # cmd4 <- paste("tabix -h",fullSS.gz,coords,">",subset_path)
  printer("TABIX:: Extracting subset of sum stats")
  dat <- data.table::fread(cmd=paste("tabix -h",fullSS.gz,coords))
  printer("++ Returning",paste(dim(dat),collapse=" x "),"data.table")
  return(dat)
}


TABIX <- function(fullSS_path,
                  subset_path,
                  is_tabix=F,
                  chrom_col="CHR",
                  position_col="POS",
                  min_POS=NA,
                  max_POS=NA,
                  chrom=NULL){
  # Check if it's already an indexed tabix file
  fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), fullSS_path, paste0(fullSS_path,".gz"))
  is_tabix <- file.exists(fullSS.gz) & file.exists(paste0(fullSS.gz,".tbi"))
  if(!is_tabix){
    fullSS.gz <- TABIX.convert_file(fullSS_path,
                                    chrom_col = chrom_col,
                                    position_col = position_col)
  } else { printer("TABIX:: Existing indexed tabix file detected") }
  # Query
  header.path <- file.path(dirname(fullSS_path),"header.txt")
  cDict <- column_dictionary(file_path = header.path)
  dat <- TABIX.query(fullSS.gz,
                      chrom=chrom,
                      start_pos=min_POS,
                      end_pos=max_POS)
  colnames(dat) <- colnames(data.table::fread(header.path))
  printer("++ Saving query ==>", subset_path)
  head(dat)
  data.table::fwrite(dat, file = subset_path, nThread = 4, sep="\t")
}




