


GITHUB.list_files <- function(creator="RajLabMSSM",
                              repo="Fine_Mapping_Shiny",
                              query=NULL,
                              return_download_api=T){
    repo_api <- file.path("https://api.github.com/repos",creator,repo,
                          "git/trees/master?recursive=1")
  req <- httr::GET(repo_api)
  httr::stop_for_status(req)
  filelist <- unlist(lapply(httr::content(req)$tree, "[", "path"), use.names = F)
  print(paste(length(filelist),"files found in GitHub repo:", file.path(creator,repo)))
  if(!is.null(query)){
    # query_string <- "*Nalls23andMe_2019.*UKB.multi_finemap.csv.gz"
    bool <- grepl(query, filelist)
    filelist <- filelist[bool]
    print(paste(length(filelist),"files found matching query."))
  }
  if(return_download_api){
    filelist <- file.path("https://github.com",creator,repo,"raw/master",filelist)
  }
  return(filelist)
}


GITHUB.download_files <- function(filelist,
                                   download_dir="./",
                                   overwrite=F,
                                   nThread=parallel::detectCores()){
  local_files <- parallel::mclapply(filelist, function(x){
    print(paste("Downloading",x))
    destfile <-  gsub("https://github.com/*.*/raw/master/www/data",
                      download_dir,x)
    dir.create(dirname(destfile), showWarnings = F, recursive = T)
    if(!file.exists(destfile) & overwrite==F) download.file(url = x, destfile=destfile)
    return(destfile)
  }, mc.cores = nThread) %>% unlist()
}


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



