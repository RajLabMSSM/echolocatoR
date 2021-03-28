


#' Install tricky R packages
#'
#' Some of \pkg{echolocatoR}'s optional R package dependencies are especially tricky to install.
#' Rather than including them in the DESCRIPTION file,
#' which would require them to install \pkg{echolocatoR} at all,
#' this functions installs them afterwards.
#' Only packages not already installed will be installed.
#'
#' Some of the main packages installed via this function include:
#' \describe{
#' \item{gaston}{CRAN}
#' \item{plotly}{CRAN}
#'
#' \item{foreign}{Archived}
#' \item{XGR}{Archived}
#' \item{refGenome}{Archived}
#'
#' \item{Rgraphviz}{Bioconductor}
#' \item{biomaRt}{Bioconductor}
#'
#' \item{knitrBootstrap}{GitHub}
#' \item{susieR}{GitHub}
#' }
#'
#' @examples
#' library(echolocatoR)
#' extra_installs()
#' @export
extra_installs <- function(cran_packages=F,
                            archived_packages=T,
                            bioc_packages=T,
                            github_packages=T){
    printer("Installing additional echolocatoR R dependencies")

    #### CRAN ####
    # current CRAN version of foreign needs R >= 4.0 - so specify legacy version
    if(cran_packages){
        message("+ Installing CRAN packages...")
        r_packages <- c("r.utils", "reticulate", "patchwork",
                        "ggrepel", "curl", "gaston", "tidyverse",
                        "crayon", "roxygen2", "coloc", "haploR", "doBy")
        # "pbmcapply", "plotly",
        required_packages <- r_packages[ ! r_packages %in% installed.packages() ]
        for(lib in required_packages){
            install.packages(lib,dependencies=TRUE)
        }
    }


    #### Archived ####
    # several packages are no longer available on CRAN - get the last approved versions
    if(archived_packages){
        message("+ Installing Archived CRAN packages...")
        if( ! "foreign" %in% installed.packages() ){
            install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz", dependencies = TRUE)
        }
        if( ! "XGR" %in% installed.packages() ){
            install.packages("https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz", dependencies = TRUE)
        }
        if( ! "refGenome" %in% installed.packages() ){
            install.packages("https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.7.tar.gz", dependencies = TRUE)
        }
        ### Above dependencies previously specified in DESCRIPTION as follows.
        ### However, removed from DESCRIPTION because these were causing serious issues with users' ability to install echolocatoR at all.
        # url::https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz,
        # url::https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.7.tar.gz,
        # url::https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz
    }

    #### Bioconductor ####
    if(bioc_packages){
        message("+ Installing Bioconductor packages...")
        library(BiocManager)
        bioc_packages <- c("supraHex", "graph", "Rgraphviz", "dnet", "rtracklayer",
                           "biomaRt", "Rsamtools", "snpStats")
        required_packages <- bioc_packages[ ! bioc_packages %in% installed.packages() ]
        for(lib in required_packages){
            BiocManager::install(lib)
        }
    }

    #### GitHub ####
    if(github_packages){
        message("+ Installing GitHub packages...")
        if(!"knitrBootstrap" %in% installed.packages() ){
            devtools::install_github('jimhester/knitrBootstrap')
        }
        if(!"susieR"  %in% installed.packages() ){
            devtools::install_github("stephenslab/susieR")
        }
    }
}

