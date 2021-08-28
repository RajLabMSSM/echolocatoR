
#### #### #### #### #### #### #### ####
#### GENERAL AND UTILITY FUNCTIONS ####
#### #### #### #### #### #### #### ####


## Documenting internal functions:
## https://www.r-bloggers.com/internal-functions-in-r-packages/

## Documenting Imports:
## https://laderast.github.io/2019/02/12/package-building-description-namespace/
# usethis::use_dev_package("tidyverse/dplyr")
# usethis::use_package("biomaRt")



#' printer
#'
#' Concatenate and print any number of items.
#'
#' @family general
#' @examples
#' n.snps <- 50
#' printer("echolocatoR::","Processing",n.snps,"SNPs...")
#' @keywords internal
#' @noRd
printer <- function(..., v=T){if(v){print(paste(...))}}



#' @family general
#' @keywords internal
.arg_list_handler <- function(arg, i){
  output <- if(length(arg)>1){arg[i]}else{arg}
  return(output)
}



#' tryCatch extension
#'
#' Extension of tryCatch function.
#'
#' @family general
#' @param input Function input.
#' @param func Function.
#' @keywords internal
tryFunc <- function(input, func, ...) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      func(input, ...)
      # The return value of `readLines()` is the actual value
      # that will be returned in case there is no condition
      # (e.g. warning or error).
      # You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", input))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", input))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message(paste("Processed URL:", input))
      message("Some other message at the end")
    }
  )
  return(out)
}







#' Replace items within DT object
#'
#' Annoyingly, there is no native function to do simple find-and-replace in the `DT` library.
#'
#' @family general
#' @return data.frame
#' @keywords internal
dt.replace <- function(DT, target, replacement){
  for(col in names(DT)) set(DT, i=which(DT[[col]]==target), j=col, value=replacement)
  return(DT)
}






# ---------------


#' @family directory functions
directory_info <- function(info_path=NULL,
                           dataset_name,
                           variable="fullSS.local"){
  Data_dirs <- data.table::fread(info_path)
  directory = subset(Data_dirs, Dataset==dataset_name, select=variable) %>% as.character()
  return(directory)
}




#' @family directory functions
#' @keywords internal
get_dataset_name <- function(file_path){
  dataset_name <- tail(strsplit(dirname(file_path), "/")[[1]],n = 1)
  return(dataset_name)
}




#' @family directory functions
#' @keywords internal
delete_subset <- function (force_new_subset, subset_path){
  # Force new file to be made
  if(force_new_subset==T){
    printer("\n + Removing existing summary stats subset...\n")
    # dataset_name <- get_dataset_name(subset_path)
    # subset_path <- paste(dirname(subset_path),"/",gene,"_",superpopulation,"_",dataset_name,"_subset.txt",sep="")
    suppressWarnings(file.remove(subset_path))
  }
}




#' Make locus-specific results folder
#'
#' @family directory functions
#' @keywords internal
#' @examples
#' locus_dir <- make_locus_dir(results_dir="./results", dataset_type="GWAS", dataset_name="Nalls23andMe_2019", locus="BST1")
make_locus_dir <- function(results_dir="./results",
                             dataset_type="dataset_type",
                             dataset_name="dataset_name",
                             locus){
  locus_dir <- file.path(results_dir, dataset_type, dataset_name, locus)
  dir.create(locus_dir, showWarnings = F, recursive = T)
  return(locus_dir)
}



#' Construct the path of the locus subset
#'
#' @family directory functions
#' @keywords internal
#' @examples
#' subset_path <- get_subset_path(results_dir="./Data/GWAS/Nalls23andMe_2019/BST1", locus="BST1")
get_subset_path <- function(subset_path="auto",
                            results_dir="./results",
                            dataset_type="dataset_type",
                            dataset_name="dataset_name",
                            locus=NULL,
                            suffix=".tsv.gz"){
  # Specify subset file name
  if(subset_path=="auto"){
    locus_dir <- make_locus_dir(results_dir = results_dir,
                                   dataset_type = dataset_type,
                                   dataset_name = dataset_name,
                                   locus = locus)
    created_sub_path <- file.path(locus_dir, paste0(locus,"_",dataset_name,"_subset",suffix) )
    return(created_sub_path)
  } else{return(subset_path)}
}



#' Extract the locus dir
#'
#' @family directory functions
#' @keywords internal
get_locus_dir <- function(subset_path){
  locus_dir <- dirname(subset_path)
  return(locus_dir)
}


#' Extract the study dir
#'
#' @family directory functions
#' @keywords internal
get_study_dir <- function(locus_dir){
  study_dir <- dirname(locus_dir)
  return(study_dir)
}



#' Create multi-finemap path
#'
#' @family directory functions
#' @keywords internal
get_multifinemap_path <- function(){
  old_file_path <- file.path(dirname(file_path),"Multi-finemap_results.txt")
}
















