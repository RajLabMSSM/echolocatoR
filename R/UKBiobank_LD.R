
#' Download LD matrices from UK Biobank
#'
#' Download pre-computed LD matrices from UK Biobank in 3Mb windows,
#' then subset to the region that overlaps with \code{subset_DT}.
#'
#' @family LD
#' @keywords internal
LD.UKBiobank <- function(subset_DT=NULL,
                         locus_dir,
                         sumstats_path=NULL,
                         chrom=NULL,
                         min_pos=NULL,

                         force_new_LD=F,
                         chimera=F,
                         server=T,
                         download_full_ld=F,
                         download_method="direct",
                         nThread=4,
                         return_matrix=F,
                         conda_env="echoR",
                         remove_tmps=T){

  LD.UKB_find_ld_prefix <- function(chrom, min_pos){
    bp_starts <- seq(1,252000001, by = 1000000)
    bp_ends <- bp_starts+3000000
    i <- max(which(bp_starts<=min_pos))
    file.name <- paste0("chr",chrom,"_", bp_starts[i],"_", bp_ends[i])
    return(file.name)
  }

  #### Support functions

  # Begin download
  if(!is.null(subset_DT)){
    finemap_dat <- subset_DT
  } else if(!is.null(sumstats_path)){
    printer("+ Assigning chrom and min_pos based on summary stats file")
    # sumstats_path="./example_data/BST1_Nalls23andMe_2019_subset.txt"
    finemap_dat <- data.table::fread(sumstats_path, nThread = 4)
  }
  chrom <- unique(finemap_dat$CHR)
  min_pos <- min(finemap_dat$POS)
  LD.prefixes <- LD.UKB_find_ld_prefix(chrom=chrom, min_pos=min_pos)
  chimera.path <- "/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD"
  alkes_url <- "https://data.broadinstitute.org/alkesgroup/UKBB_LD"
  URL <- alkes_url

  printer("+ UKB LD file name:",LD.prefixes)
  LD_dir <- file.path(locus_dir, "LD")
  dir.create(LD_dir, showWarnings = F, recursive = T)
  locus <- basename(locus_dir)
  UKBB.LD.RDS <- file.path(LD_dir, paste0(locus,".UKB_LD.RDS"))

  if(file.exists(UKBB.LD.RDS) & force_new_LD==F){
    printer("+ LD:: Pre-existing UKB_LD.RDS file detected. Importing...")
    LD_matrix <- readRDS(UKBB.LD.RDS)
  } else {
    if(download_method!="direct"){
      if(download_full_ld | force_new_LD | download_method %in% c("wget","axel")){
        printer("+ LD:: Downloading full .gz/.npz UKB files and saving to disk.")
        URL <- LD.download_UKB_LD(LD.prefixes = LD.prefixes,
                                  locus_dir = locus_dir,
                                  background = F,
                                  force_overwrite = force_new_LD,
                                  download_method = download_method)
        server <- F
      } else {
        if(chimera){
          if(file.exists(file.path(chimera.path, paste0(LD.prefixes,".gz")))  &
             file.exists(file.path(chimera.path, paste0(LD.prefixes,".npz"))) ){
            printer("+ LD:: Pre-existing UKB LD gz/npz files detected. Importing...")
            URL <- chimera.path
          }
        } else {
          URL <- file.path(URL, LD.prefixes)
          printer("+ LD:: Importing UKB LD file directly to R from:")
          print(URL)
        }
      }
    } else {
      URL <- file.path(URL, LD.prefixes)
      printer("+ LD:: Importing UKB LD file directly to R from:")
      print(URL)
    }



    # RSIDs file
    # rsids <- data.table::fread(gz.path, nThread = 4)
    printer("+ LD:: ...this could take some time...")
    CONDA.activate_env(conda_env = conda_env)
    reticulate::source_python(system.file("tools","load_ld.py",package = "echolocatoR"))

    # load_ld(ld_prefix=URL, server=F)
    ld.out <- tryFunc(input = URL, load_ld, server)
    # LD matrix
    ld_R <- ld.out[[1]]
    # head(ld_R)[1:10,]
    # SNP info
    # ld_snps <- data.table::data.table( reticulate::py_to_r(ld.out[[2]]) )
    ld_snps <- data.table::data.table(ld.out[[2]])
    row.names(ld_R) <- ld_snps$rsid
    colnames(ld_R) <- ld_snps$rsid

    # remove(ld.out)
    # ld_snps.sub <- subset(ld_snps, position %in% finemap_dat$POS)
    indices <- which(ld_snps$position %in% finemap_dat$POS)
    ld_snps.sub <- ld_snps[indices,]
    LD_matrix <- ld_R[indices, indices]
    row.names(LD_matrix) <- ld_snps.sub$rsid
    colnames(LD_matrix) <- ld_snps.sub$rsid
    LD_matrix[is.na(LD_matrix)] <- 0
    # Save LD matrix as RDS
    printer("LD matrix dimensions", paste(dim(LD_matrix),collapse=" x "))
    printer("+ LD:: Saving LD =>",UKBB.LD.RDS)
    dir.create(dirname(UKBB.LD.RDS), showWarnings = F, recursive = T)
    saveRDS(LD_matrix, UKBB.LD.RDS)

    if(remove_tmps){
      printer("+ Removing .gz/.npz files.")
      if(file.exists(paste0(URL,".gz"))){ file.remove(paste0(URL,".gz")) }
      if(file.exists(paste0(URL,".npz"))){ file.remove(paste0(URL,".npz")) }
    }
  }
  if(return_matrix){
    return(LD_matrix)
  } else {
    return(UKBB.LD.RDS)
  }
}




LD.download_UKB_LD <- function(LD.prefixes,
                               locus_dir,
                               alkes_url="https://data.broadinstitute.org/alkesgroup/UKBB_LD",
                               background=T,
                               force_overwrite=F,
                               download_method="direct",
                               nThread=4){
  for(f in LD.prefixes){
    gz.url <- file.path(alkes_url,paste0(f,".gz"))
    npz.url <- file.path(alkes_url,paste0(f,".npz"))

    for(furl in c(gz.url, npz.url)){
      if(tolower(download_method)=="axel"){
        out.file <- axel(input_url = furl,
                         output_path = file.path(locus_dir,"LD"),
                         background = background,
                         nThread = nThread,
                         force_overwrite = force_overwrite)
      }
      if(tolower(download_method)=="wget"){
        out.file <- wget(input_url = furl,
                         output_path = file.path(locus_dir,"LD"),
                         background = background,
                         force_overwrite = force_overwrite)
      }
    }
  }
  return(gsub("*.npz$","",out.file))
}



#' Convert .RDS file back to .npz format
#'
#' @family LD
#' @keywords internal
#' @examples
#' \dontrun{
#' data("BST1")
#' npz_path <- LD.rds_to_npz(rds_path="/Users/schilder/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/BST1/plink/UKB_LD.RDS")
#' }
LD.rds_to_npz <- function(rds_path,
                          conda_env="echoR",
                          verbose=T){
  printer("POLYFUN:: Converting LD .RDS to .npz:",rds_file, v=verbose)
  LD_matrix <- readRDS(rds_path)
  reticulate::use_condaenv(condaenv = conda_env)
  np <- reticulate::import(module = "numpy")
  npz_path <- gsub(".RDS",".npz", rds_path)
  np$savez(npz_path, as.matrix(LD_matrix))
  return(npz_path)
}

