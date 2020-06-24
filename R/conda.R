
# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~CONDA~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~# ~~~~~~~~~~~~~




#' Install conda if it's missing
#'
#' @family conda
CONDA.install <- function(conda_path="auto"){
  conda_version <- NULL
  try({conda_version <- reticulate::conda_version(conda = conda_path)})
  if(is.null(conda_version)){
    printer("+ CONDA:: conda not detected. Installing with reticulate...")
    reticulate::conda_install()
  }
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



#' Identify the right conda env to use
#'
#' @family conda
CONDA.find_env_path <- function(conda_env="echoR"){
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

#' Install necessary command line and python tools
#'
#' \describe{
#' \item{plink}{https://anaconda.org/bioconda/plink}
#' \item{tabix}{https://anaconda.org/bioconda/tabix}
#' }
#'
#' @family conda
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
  CONDA.activate_env(conda_env=envname)
}






#' Create conda env from yaml file
#' @keywords internal
#' @family conda
CONDA.env_from_yaml <- function(yaml_path=system.file("conda","echoR.yml",package = "echolocatoR")){
  cmd <- paste("conda env create -f",yaml_path)
  print(cmd)
  system(cmd)
}




#' Create new conda env from list of dependencies
#' @keywords internal
#' @family conda
CONDA.env_from_list <- function(libraries = c("numpy",
                                                  "scipy",
                                                  "scikit-learn",
                                                  "pandas",
                                                  "tqdm",
                                                  "pyarrow",
                                                  "bitarray",
                                                  "networkx",
                                                  "rpy2",
                                                  "r-ckmeans.1d.dp")){
  # NOTE: version specification must use quotes
  cmd <- paste("conda create -n polyfun_venv python=3.7.3",
               "numpy scipy scikit-learn 'pandas>=0.25.0'",
               paste(libraries, collapse=" "))
  print(cmd)
}




