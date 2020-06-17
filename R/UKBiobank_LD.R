
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
                         nThreads=4,
                         return_matrix=F,
                         remove_tmps=T){

  LD.UKB_find_ld_prefix <- function(chrom, min_pos){
    bp_starts <- seq(1,252000001, by = 1000000)
    bp_ends <- bp_starts+3000000
    i <- max(which(bp_starts<=min_pos))
    file.name <- paste0("chr",chrom,"_", bp_starts[i],"_", bp_ends[i])
    return(file.name)
  }


  LD.download_UKB_LD <- function(LD.prefixes,
                                 locus_dir,
                                 alkes_url="https://data.broadinstitute.org/alkesgroup/UKBB_LD",
                                 background=T,
                                 force_overwrite=F,
                                 download_method="direct",
                                 nThreads=4){
    for(f in LD.prefixes){
      gz.url <- file.path(alkes_url,paste0(f,".gz"))
      npz.url <- file.path(alkes_url,paste0(f,".npz"))

      for(furl in c(gz.url, npz.url)){
        if(tolower(download_method)=="axel"){
          out.file <- axel(input_url = furl,
                           output_path = file.path(locus_dir,"LD"),
                           background = background,
                           nThreads = nThreads,
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





# LD.UKBiobank_multi_download <- function(out.path = "/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD/"){
#   # Download all UKBB LD files
#   # wget -r -np -A '*.gz' https://data.broadinstitute.org/alkesgroup/UKBB_LD/
#   # wget -r -np -A '*.npz' https://data.broadinstitute.org/alkesgroup/UKBB_LD/
#
#   # Figure out the names of LD files you'll need
#   merged_DT <- merge_finemapping_results(minimum_support = 0)
#   locus_coords <- merged_DT %>% dplyr::group_by(Gene) %>%
#     dplyr::summarise(chrom=unique(CHR), min_pos=min(POS), max_pos=max(POS))
#   LD.file.list <- lapply(1:nrow(locus_coords), function(i){
#     return(LD.find_ld_prefix(chrom = locus_coords$chrom[i], min_pos=locus_coords$min_pos[i]))
#   }) %>% unlist()
#
#   LD.download_UKB_LD(LD.file.list = LD.file.list,
#                      out.path = out.path)
#   return(LD.file.list)
# }


# POLYFUN.load_conda <- function(server=F){
#   printer("POLYFUN:: Activating polyfun_venv...")
#   if(server){
#     reticulate::use_condaenv("polyfun_venv")
#   } else {
#     reticulate::use_condaenv("polyfun_venv", conda = "/usr/local/anaconda3/condabin/conda")
#     # conda_list("/usr/local/anaconda3/condabin/conda")
#   }
# }
#



# LD.UKB_find_ld_prefix <- function(chrom, min_pos){
#   bp_starts <- seq(1,252000001, by = 1000000)
#   bp_ends <- bp_starts+3000000
#   i <- max(which(bp_starts<=min_pos))
#   file.name <- paste0("chr",chrom,"_", bp_starts[i],"_", bp_ends[i])
#   return(file.name)
# }
#
#
#
#
# LD.UKBiobank <- function(finemap_dat,
#                          subset_path,
#                          force_new_LD=F,
#                          polyfun="./echolocatoR/tools/polyfun",
#                          server=F){
#   alkes_url <- "https://data.broadinstitute.org/alkesgroup/UKBB_LD"
#   URL <- alkes_url
#   chrom <- unique(finemap_dat$CHR)
#   min_pos <- min(finemap_dat$POS)
#   file.name <- LD.UKB_find_ld_prefix(chrom=chrom, min_pos=min_pos)
#   printer("+ UKB LD file name:",file.name)
#
#   gz.path <- file.path(alkes_url, paste0(file.name,".gz"))
#   npz.path <- file.path(alkes_url, paste0(file.name,".npz"))
#
#   chimera.path <- file.path("/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD")
#   UKBB.LD.file <- file.path(subset_path,"plink/UKB_LD.RDS")
#
#   if(server){
#     if(file.exists(file.path(chimera.path, paste0(file.name,".gz")))  &
#        file.exists(file.path(chimera.path, paste0(file.name,".npz"))) ){
#       printer("+ LD:: Pre-existing UKB LD gz/npz files detected. Importing...")
#       URL <- chimera.path
#     } else {
#       printer("+ LD:: Downloading .gz/.npz and saving to disk.")
#       LD.download_UKB_LD(LD.file.list = file.name,
#                          out.path = URL,
#                          background = F)
#     }
#   } else {
#     printer("+ LD:: Importing UKB LD file from alkesgroup database directly to R.")
#   }
#
#   # RSIDs file
#   # rsids <- data.table::fread(gz.path, nThread = 4)
#   if(file.exists(UKBB.LD.file) & force_new_LD!=T){
#     printer("POLYFUN:: Pre-existing UKB_LD.RDS file detected. Importing...")
#     LD_matrix <- readRDS(UKBB.LD.file)
#   } else {
#     POLYFUN.load_conda(server = server)
#     reticulate::source_python(file.path(polyfun,"load_ld.py"))
#     printer("+ LD:: ...this could take some time...")
#
#     ld.out <- tryFunc(input = file.path(URL, file.name), load_ld)
#     if(is.na(ld.out)){
#       ld.out <- load_ld(ld_prefix = file.path(URL, file.name), npz_suffix='2')
#     }
#     # LD matrix
#     ld_R <- ld.out[[1]]
#     # head(ld_R)[1:10]
#     # SNP info
#     # ld_snps <- data.table::data.table( reticulate::py_to_r(ld.out[[2]]) )
#     ld_snps <- data.table::data.table(ld.out[[2]])
#     row.names(ld_R) <- ld_snps$rsid
#     colnames(ld_R) <- ld_snps$rsid
#
#     # remove(ld.out)
#     # ld_snps.sub <- subset(ld_snps, position %in% finemap_dat$POS)
#     indices <- which(ld_snps$position %in% finemap_dat$POS)
#     ld_snps.sub <- ld_snps[indices,]
#     LD_matrix <- ld_R[indices, indices]
#     row.names(LD_matrix) <- ld_snps.sub$rsid
#     colnames(LD_matrix) <- ld_snps.sub$rsid
#     LD_matrix[is.na(LD_matrix)] <- 0
#     printer("LD matrix dimensions", paste(dim(LD_matrix),collapse=" x "))
#     printer("+ POLYFUN:: Saving LD =>",UKBB.LD.file)
#     dir.create(dirname(UKBB.LD.file), showWarnings = F, recursive = T)
#     saveRDS(LD_matrix, UKBB.LD.file)
#   }
#   return(LD_matrix)
# }
#
#
# LD.list_all_LDfiles <- function(alkes_ld_excel="./echolocatoR/tools/Alkes_UKB_LD.xlsx"){
#   alkes_url="https://data.broadinstitute.org/alkesgroup/UKBB_LD"
#   xl <- readxl::read_excel(alkes_ld_excel)
#   npz.paths <- file.path(alkes_url, grep(".npz$",xl$Name, value = T))
#   paste(paste0("'",npz.paths,"'"), collapse=" ")
#   return(npz.paths)
# }
#
# LD.download_UKB_LD <- function(LD.file.list,
#                                out.path="/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD/",
#                                alkes_url="https://data.broadinstitute.org/alkesgroup/UKBB_LD",
#                                background=T,
#                                method="wget"){
#   flags <- ifelse(background,"-bqc","qc")
#   for(f in LD.file.list){
#     gz.path <- file.path(alkes_url,paste0(f,".gz"))
#     npz.path <- file.path(alkes_url,paste0(f,".npz"))
#     # https://stackoverflow.com/questions/21365251/how-to-run-wget-in-background-for-an-unattended-download-of-files
#     if(method=="wget"){
#       ## -bqc makes wget run in the background quietly
#       # gz file
#       system(paste("wget",gz.path,"-np",flags,"-P",out.path))
#       # npz file
#       cmd <- paste("wget",npz.path,"-np",flags,"-P",out.path)
#       print(cmd)
#       system(cmd)
#     } else if (method=="aria2") {
#       printer("++ Enabling parallelized download with aria2")
#       system(paste("aria2c -d",out.path, gz.path))
#       cmd <- paste("aria2c -d",out.path, gz.path)
#       print(cmd)
#       system(cmd)
#     }
#   }
# }
#
# # LD.UKBiobank_multi_download <- function(out.path = "/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD/"){
# #   # Download all UKBB LD files
# #   # wget -r -np -A '*.gz' https://data.broadinstitute.org/alkesgroup/UKBB_LD/
# #   # wget -r -np -A '*.npz' https://data.broadinstitute.org/alkesgroup/UKBB_LD/
# #
# #   # Figure out the names of LD files you'll need
# #   merged_DT <- merge_finemapping_results(minimum_support = 0)
# #   locus_coords <- merged_DT %>% dplyr::group_by(Gene) %>%
# #     dplyr::summarise(chrom=unique(CHR), min_pos=min(POS), max_pos=max(POS))
# #   LD.file.list <- lapply(1:nrow(locus_coords), function(i){
# #     return(LD.find_ld_prefix(chrom = locus_coords$chrom[i], min_pos=locus_coords$min_pos[i]))
# #   }) %>% unlist()
# #
# #   LD.download_UKB_LD(LD.file.list = LD.file.list,
# #                      out.path = out.path)
# #   return(LD.file.list)
# # }
#
#
#
# POLYFUN.load_conda <- function(server=F){
#   printer("POLYFUN:: Activating polyfun_venv...")
#   if(server){
#     reticulate::use_condaenv("polyfun_venv")
#   } else {
#     reticulate::use_condaenv("polyfun_venv", conda = "/usr/local/anaconda3/condabin/conda")
#     # conda_list("/usr/local/anaconda3/condabin/conda")
#   }
# }
#
#
