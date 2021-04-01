


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
#' \item{XGR}{xgr}
#' \item{foreign}{xgr}
#' \item{refGenome}{xgr}
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
                           xgr=T,
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
        required_packages <- r_packages[ ! r_packages %in% row.names(installed.packages()) ]
        for(lib in required_packages){
            install.packages(lib, dependencies=T)
        }
    }


    #### XGR ####
    ## XGR has been extra tricky to install and
    ## seems to be very sensitive to R version and the versions of its dependencies.
    ## So employ a tiered approach to installation: Bioconductor ==> CRAN ==> archived
    # several packages are no longer available on CRAN - get the last approved versions
    if(xgr){
        #### Try via Bioconductor ####
        message("+ Attempt 1: Installing XGR and deps via Bioconductor...")
        if( ! "XGR" %in% row.names(installed.packages()) ){
            try({
                BiocManager::install("hfang-bristol/XGR", dependencies=T)
            })
        }
        #### Attempt 2: Try via CRAN ####
        message("+ Attempt 2: Installing XGR and deps via CRAN...")
        if( ! "XGR" %in% row.names(installed.packages()) ){
            try({
                install.packages("XGR", dependencies = T)
            })
        }
        #### Attempt 3: Try via archived ####
        if( ! "XGR" %in% row.names(installed.packages()) ){
            #### Archived version ####
            message("+ Attempt 3: Installing XGR and deps via Archives...")
            if( ! "foreign" %in% row.names(installed.packages()) ){
                install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz", dependencies = T)
            }
            if( ! "refGenome" %in% row.names(installed.packages()) ){
                install.packages("https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.7.tar.gz", dependencies = T)
            }

            install.packages("https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz", dependencies = T)
            ### Above dependencies previously specified in DESCRIPTION as follows.
            ### However, removed from DESCRIPTION because these were causing serious issues with users' ability to install echolocatoR at all.
            # url::https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz,
            # url::https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.7.tar.gz,
            # url::https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz
        }
    }

    #### Bioconductor ####
    if(bioc_packages){
        message("+ Installing Bioconductor packages...")
        library(BiocManager)
        bioc_packages <- c("supraHex", "graph", "Rgraphviz", "dnet", "rtracklayer",
                           "biomaRt", "Rsamtools", "snpStats")
        required_packages <- bioc_packages[ ! bioc_packages %in% row.names(installed.packages()) ]
        for(lib in required_packages){
            BiocManager::install(lib)
        }
    }

    #### GitHub ####
    if(github_packages){
        message("+ Installing GitHub packages...")
        if(!"knitrBootstrap" %in% row.names(installed.packages()) ){
            devtools::install_github('jimhester/knitrBootstrap')
        }
        if(!"susieR"  %in% row.names(installed.packages() )){
            devtools::install_github("stephenslab/susieR")
        }
    }
}

