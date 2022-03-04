
#' echolocatoR Fine-mapping portal: metadata
GITHUB.portal_metadata <- function(verbose=T){
  printer("Fetching echolocatoR Fine-mapping Portal study metadata.",v=verbose)
  meta <- data.table::fread("https://github.com/RajLabMSSM/Fine_Mapping_Shiny/raw/master/www/metadata/study_metadata.csv.gz")
  return(meta)
}





#' Search and download fine-mapping files
#'
#' Search the \href{https://github.com/RajLabMSSM/Fine_Mapping_Shiny}{echolocatoR Fine-mapping Portal}
#' for fine-mapping results, LD, and locus plots.
#' @export
#' @examples
#' local_finemap <- GITHUB.portal_query(dataset_types="GWAS",
#'                                      phenotypes = c("schizophrenia","parkinson"),
#'                                      file_types = "multi_finemap",
#'                                      loci = c("BST1","CHRNB1","LRRK2",1:3),
#'                                      LD_panels=c("UKB","1KGphase3"))
GITHUB.portal_query <- function(dataset_types=NULL,
                                datasets=NULL,
                                phenotypes=NULL,
                                loci=NULL,
                                LD_panels=c("UKB","1KGphase1","1KGphase3"),
                                file_types=c("multi_finemap","LD","plot"),
                                results_dir=tempdir(),
                                overwrite=F,
                                nThread=parallel::detectCores()-2,
                                verbose=T){
  #### Search metadata ####
  meta <- GITHUB.portal_metadata(verbose=verbose)
  meta <- if(!is.null(dataset_types)) subset(meta, tolower(dataset_type) %in% tolower(dataset_types)) else meta
  meta <- if(!is.null(datasets)) meta[grepl(paste(datasets,collapse = "|"),meta$dataset,ignore.case = T),] else meta
  meta <- if(!is.null(phenotypes)) meta[grepl(paste(phenotypes,collapse = "|"),meta$phenotype,ignore.case = T),] else meta
  printer("+",nrow(meta),"datasets remain after filtering.",v=verbose)

  #### Find URLs ####
  file_type_dict <- c("multi_finemap"=".multi_finemap.csv.gz",
                      "LD"=".LD.csv.gz",
                      "plot"=".png")
  file_urls <- lapply(file_types, function(ftype){
   printer("+ Searching for",ftype,"files...",v=verbose)
   remote_finemap <- GITHUB.list_files(creator="RajLabMSSM",
                                       repo="Fine_Mapping_Shiny",
                                       query=file_type_dict[[ftype]],
                                       branch = "master", #IMPORTANT! not "main" like other repos for some reason
                                       # query= paste(paste0("(",paste(unique(meta$dataset_type), collapse = "|"),")"),
                                       #              paste0("(",file_type_dict[[ftype]],")"),sep="&&"),
                                       verbose = F)
   return(data.table::data.table(URL=remote_finemap, file_type=ftype))
 }) %>% data.table::rbindlist() %>%
   # make sure to remove anything that's not in the data folder (e.g. icons)
   subset(startsWith(URL,"https://github.com/RajLabMSSM/Fine_Mapping_Shiny/raw/master/www/data"))

  #### Format url results ####
  file_filt <- file_urls %>%
   dplyr::mutate(url_strip=gsub("https://github.com/RajLabMSSM/Fine_Mapping_Shiny/raw/master/www/data/","",URL)) %>%
   tidyr::separate(col = "url_strip",sep = "/", into = c("dataset_type","dataset","locus"), remove = F, extra = "drop") %>%
   subset(dataset_type %in% unique(meta$dataset_type) &
          dataset %in% unique(meta$dataset))

 ## Let's just download the loci of interest.
 if(!is.null(loci)) file_filt <- subset(file_filt, locus %in% loci)
 printer("+",nrow(file_filt),"unique files identified.",v=verbose)
 ## Filter by LD panel
 file_filt <- file_filt[grepl(paste(LD_panels,collapse = "|"),file_filt$URL),]

  ### Download files ####
  local_finemap <- GITHUB.download_files(filelist = unique(file_filt$URL),
                                         download_dir = results_dir,
                                         overwrite = overwrite,
                                         nThread = nThread,
                                         verbose = verbose)
 printer("+ Returning local file paths.",v=verbose)
 return(local_finemap)
}




GITHUB.list_files <- function(creator="neurogenomics",
                              repo="MAGMA_Files",
                              branch=c("main","master"),
                              query=NULL,
                              return_download_api=T,
                              verbose=T){
  repo_api <- file.path("https://api.github.com/repos",creator,repo,
                        paste0("git/trees/",branch[1],"?recursive=1"))
  req <- httr::GET(repo_api)
  httr::stop_for_status(req)
  filelist <- unlist(lapply(httr::content(req)$tree, "[", "path"), use.names = F)
  printer(paste(length(filelist),"files found in GitHub repo:", file.path(creator,repo)),v=verbose)
  if(!is.null(query)){
    # query_string <- "*Nalls23andMe_2019.*UKB.multi_finemap.csv.gz"
    bool <- grepl(query, filelist)
    filelist <- filelist[bool]
    printer(paste(length(filelist),"files found matching query."),v=verbose)
  }
  if(return_download_api){
    filelist <- file.path("https://github.com",creator,repo,"raw",branch,filelist)
  }
  return(filelist)
}


GITHUB.download_files <- function(filelist,
                                  download_dir=tempdir(),
                                  overwrite=F,
                                  nThread=parallel::detectCores()-2,
                                  verbose=T){
  printer("+ Downloading",length(filelist),"files...",v=verbose)
  local_files <- unlist(parallel::mclapply(filelist, function(x){
    print(paste("Downloading",x))
    branch <- stringr::str_split(string = x, pattern = "/")[[1]][7]
    folder_structure <- paste(stringr::str_split(string = x, pattern = "/")[[1]][-c(1:7)], collapse="/")
    destfile <- gsub("/www/data","",file.path(download_dir, folder_structure))
    dir.create(dirname(destfile), showWarnings = F, recursive = T)
    if(!file.exists(destfile) & overwrite==F) download.file(url = x, destfile=destfile)
    return(destfile)
  }, mc.cores = nThread))
  return(local_files)
}


#' echolocatoR Fine-mapping portal: data dictionary
#'
#' @export
GITHUB.make_data_dict <- function(named_lists){
  data_dict <- list()
  for(x in names(named_list)){
    dat <- named_list[[x]]
    data_dict[[x]] <- as.list(setNames(dat, basename(dirname(dirname(dat))) ) )
  }
  # Add the locus dir as a bonus
  data_dict$locus_dir <- as.list(setNames(dirname(dirname(dat)), basename(dirname(dirname(dat))) ))
  # Dataset dir
  data_dict$dataset_dir <- dirname(dirname(dirname(dat)))[1]
  # Dataset name
  data_dict$dataset_type <- basename(dirname(dirname(dirname(dirname(dat)))))[1]
  # data
  data_dict$dataset <- basename(dirname(dirname(dirname(dat))))[1]
  return(data_dict)
}




GITHUB.find_pages <- function(creator="RajLabMSSM",
                              repo="Fine_Mapping",
                              local_repo=NULL,
                              return_table=T,
                              save_path=NULL){
  if(is.null(local_repo)){
    filelist <- GITHUB.list_files(creator=creator,
                                  repo=repo,
                                  query="*.*.html",
                                  return_download_api=F)
  } else {
    filelist <- gsub("^[.][/]","",list.files(path = local_repo, pattern="*.*.html", full.names=T, recursive=T))
  }
  gh_pages_url <- file.path(paste0("https://",creator,".github.io"),repo)
  gh_pages_links <- file.path(gh_pages_url, filelist)
  if(return_table){
    links_df <- data.frame(creator=creator,
                           repo=repo,
                           dir=unlist(lapply(filelist, function(x){strsplit(x,"/")[[1]][1]} )),
                           subdir=unlist(lapply(filelist, function(x){strsplit(x,"/")[[1]][2]} )),
                           url=gh_pages_links,
                           link=paste0("<a href='",gh_pages_links,"' target='blank'>",filelist,"</a>"))
    if(!is.null(save_path)){
      print(paste("Writing links to ==>",save_path))
      data.table::fwrite(links_df,save_path, sep=",")
      write.table(paste(pages[["link"]],collapse="\n\n"), "results.md")
    }
    return(links_df)
  } else {
    return(gh_pages_links)
  }
}



