#' source_all
#'
#' Source all files in a directory at once.
#' Also loads selected libraries.
#'
#' @inheritParams base::list.files
#' @keywords internal
source_all <- function(path="R/",
                       pattern="*.R$",
                       packages="dplyr"){
    for(x in packages){
        library(x, character.only=TRUE)
    }
    ### Source all internal funcs at once
    file.sources = list.files(path =path,
                              pattern = pattern,
                              full.names = TRUE, ignore.case = TRUE)
    message("Sourcing ",length(file.sources)," files.")
    out <- sapply(file.sources,source)
}
