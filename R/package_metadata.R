package_metadata <- function(pkgs = "echoverse",
                             fields = c("Package",
                                         "Title",
                                         "Description",
                                         "URL",
                                         "Version",
                                         "Depends",
                                         "Imports",
                                         "Suggests",
                                         "Remotes",
                                         "SystemRequirements"),
                              verbose = FALSE){
  echoverse <- c("echolocatoR",
                 "echodata",
                 "echotabix",
                 "echoannot",
                 "echoconda",
                 "echoLD",
                 "echoplot",
                 "catalogueR",
                 "echofinemap",
                 "downloadR")
  if(tolower(pkgs)=="echoverse"){
    messager("Collecting metadata for all echoverse modules.",v=verbose)
    pkgs <- echoverse
  } else {
    messager("Collecting metadata for",length(pkgs),"packages",v=verbose)
  }
  #### Split func ####
  parse_deps <- function(d,
                         field,
                         split=","){
    gsub("\n","",strsplit(d[[field]],split=split)[[1]])
  }
  #### Iterate ####
  meta <- lapply(pkgs, function(pkg){
    tryCatch({
      messager(pkg, v=verbose)
      d <- utils::packageDescription(pkg)
      data.table::data.table(
        t(
          lapply(fields,function(x){
            messager(" -- ",x,v=verbose)
            if(is.null(d[[x]]))  return(NULL)
            parse_deps(
              d = d,
              field = x,
              split = if(x %in% c("Title","Description")) "______" else ",")
          }) %>% `names<-`(fields)
        )
      )
    }, error = function(e){warning(e); NULL})
  }) %>% data.table::rbindlist()
  meta$Package <- unlist(meta$Package)
  data.table::setkey(meta,Package)
  return(meta)
}
