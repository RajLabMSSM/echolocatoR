
# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~CONDA~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~



#' @section conda functions:
#' Functions for setting up \emph{conda} environments with R, Python, and command line tools.



#' Install conda if it's missing
#'
#' @family  conda functions
CONDA.install <- function(){
  conda_version <- NULL
  try({conda_version <- reticulate::conda_version()})
  if(is.null(conda_version)){
    reticulate::conda_install()
  }
}

#' Activate conda env
#'
#' @family conda functions
CONDA.activate_env <- function(condaenv="echoR"){
  reticulate::use_condaenv(condaenv = condaenv)
}


#' Install necessary command line and python tools
#'
#' \describe{
#' \item{plink}{https://anaconda.org/bioconda/plink}
#' \item{tabix}{https://anaconda.org/bioconda/tabix}
#' }
#'
#' @family conda functions
#' @param envname The conda environment where you want to install \emph{plink}.
#' By default uses \emph{echoR}, the conda environment distributed with \emph{echolocatoR}.
CONDA.create_env <- function(envname="echoR",
                                   packages=c("plink","tabix",
                                              # python,
                                              "pandas")){
  envs <- reticulate::conda_list()$name
  if(!envname %in% envs){
    reticulate::conda_install(envname = envname,
                              packages = packages,
                              channel = "bioconda",
                              python_version = 3.7)
  }
  CONDA.activate_env(condaenv=envname)
}


