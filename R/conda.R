
# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~CONDA~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~


printer <- function(..., v=T){if(v){print(paste(...))}}



#' Install reticulate
#'
#' \emph{reticulate} often doesn't install very well via CRAN.
#' This function helps do it correctly.
CONDA.install_reticulate <- function(dependencies=c("devtools",
                                                    "reticulate",
                                                    "lattice",
                                                    "jsonlite",
                                                    "Matrix",
                                                    "rappdirs",
                                                    "Rcpp")){
  if("reticulate" %in% installed.packages()){
    print("CONDA:: `reticulate` already installed.")
  } else{
    required_packages <- dependencies[ ! dependencies %in% installed.packages() ]
    if(length(required_packages)>0){
      for(lib in required_packages){
        install.packages(lib, dependencies=T)
      }
    }
  }
}







#' Install conda if it's missing
#'
#' @family conda
CONDA.install <- function(conda_path="auto"){
  conda_version <- NULL
  try({conda_version <- reticulate::conda_version(conda = conda_path)})
  if(is.null(conda_version)){
    printer("+ CONDA:: conda not detected. Installing with reticulate...")
    reticulate::conda_install()
  } else {printer("+ CONDA:: conda already installed.")}
}




#' Activate conda env
#'
#' @family conda
#' @examples
#' CONDA.activate_env(conda_env="echoR")
CONDA.activate_env <- function(conda_env="echoR"){
  CONDA.install()
  env_list <- reticulate::conda_list()
  if(conda_env %in% env_list$name){
    printer("+ CONDA:: Activating conda env",paste0("'",conda_env,"'"))
    reticulate::use_condaenv(condaenv = conda_env)
  } else {
    printer("+ CONDA::",paste0("'",conda_env,"'"),"conda environment not found. Using default 'base' instead.")
    reticulate::use_condaenv(condaenv = "base")
  }
}



#' Find the python file for a specific env
#'
#' @family conda
CONDA.find_python_path <- function(conda_env="echoR"){
  CONDA.install()
  env_list <- reticulate::conda_list()
  if(conda_env %in% env_list$name){
    python_path <- subset(env_list,  name==conda_env)$python
  } else {
    printer("+ CONDA::",paste0("'",conda_env,"'"),"conda environment not found. Using default 'python' instead.")
    python_path <- "python"
  }
 return(python_path)
}




#' Find the R library for a specific env
#'
#' @family conda
CONDA.find_env_Rlib <- function(conda_env="echoR"){
 conda_path <- dirname(dirname(CONDA.find_python_path(conda_env = conda_env)))
 env_Rlib <- file.path(conda_path,"lib/R/library/")
 return(env_Rlib)
}




#' Create conda env for \emph{echolocatoR}
#'
#' Create a new env (or update and existing one)
#' with the necessary Python, R, and command line packages
#' to run \emph{echolocatoR}.
#'
#' \describe{
#' \item{plink}{https://anaconda.org/bioconda/plink}
#' \item{tabix}{https://anaconda.org/bioconda/tabix}
#' }
#' @family conda
#' @param envname The conda environment where you want to install \emph{plink}.
#' By default uses \emph{echoR}, the conda environment distributed with \emph{echolocatoR}.
CONDA.create_echoR_env <- function(conda_env="echoR",
                                   python_version=NULL,
                                   channels=c("conda-forge","bioconda","r"),
                                   python_packages=c("pandas>=0.25.0",
                                                     "pyarrow",
                                                     "scikit-learn",
                                                     "bitarray",
                                                     "networkx",
                                                     "rpy2",
                                                     "scipy",
                                                     "pandas-plink"),
                                   r_packages=c("r-base",
                                                # "r>=3.6.3",
                                                "r-devtools"
                                                # "r-biocmanager",
                                                # "r-reticulate",
                                                # "r-data.table",
                                                # "r-ggplot2",
                                                # "r-wavethresh",
                                                # "r-lattice",
                                                # "r-ckmeans.1d.dp",
                                                # "r-stringi",
                                                # "r-matrixstats",
                                                # "r-expm",
                                                # "r-rlang",
                                                # "r-xgr",
                                                ),
                                   cli_packages=c("tabix",
                                                  "plink",
                                                  "macs2"),
                                   force_install=F,
                                   auth_token=devtools::github_pat()){
  # Make sure conda is installed to begin with
  CONDA.install()
  # conda_path <- reticulate::conda_binary()
  envs <- reticulate::conda_list()$name
  if(!envname %in% envs | force_install){
    reticulate::conda_install(envname = conda_env,
                              packages = c(python_packages, r_packages, cli_packages),
                              channel = channels,
                              python_version = python_version)
  }
  CONDA.activate_env(conda_env=conda_env)
  # Get the path to the conda env
  env_path <- dirname(dirname(CONDA.find_python_path(conda_env=conda_env)))
  env_Rlib <- CONDA.find_env_Rlib(conda_env = conda_env)
  # Have echolocatoR install itself into the new env
  devtools::install_github(repo="RajLabMSSM/echolocatoR",
                           auth_token=auth_token,
                           upgrade="always",
                           lib=env_Rlib, force = T)

  return(env_path)
}






#' Create conda env from yaml file
#' @keywords internal
#' @family conda
CONDA.env_from_yaml <- function(yaml_path=system.file("conda","echoR.yml",package = "echolocatoR")){
  cmd <- paste("conda env create -f",yaml_path)
  print(cmd)
  system(cmd)
}






