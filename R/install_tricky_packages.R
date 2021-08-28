#' Install tricky packages
#'
#' Some packages don't install very well via the DESCRIPTION file
#' (e.g. wrong versions, wrong sources).
#' This function ensures they're actually installed properly.
#' @export
install_tricky_packages <- function(){
  if(!"foreign" %in% installed.packages()){
    install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz",
                     dependencies = T)
  }
  if(!"XGR" %in% installed.packages()){
    install.packages("https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz",
                     dependencies = T)
  }
  if(!"refGenome" %in% installed.packages()){
    install.packages("https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.7.tar.gz",
                     dependencies = T)
  }
}

