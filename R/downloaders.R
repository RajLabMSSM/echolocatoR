


#' wget
#'
#' R wrapper for wget
#' @family download functions
#' @keywords internal
wget <- function(input_url,
                 output_path,
                 background=T,
                 force_overwrite=F,
                 quiet=F,
                 show_progress=T,
                 continue=T){
  # https://stackoverflow.com/questions/21365251/how-to-run-wget-in-background-for-an-unattended-download-of-files
  ## -bqc makes wget run in the background quietly
  dir.create(output_path, showWarnings = F, recursive = T)
  out.file <- file.path(output_path,basename(input_url))
  cmd <- paste("wget",input_url,
               "-np",
               ifelse(background,"-b",""),
               ifelse(continue,"-c",""),
               ifelse(quiet,"-q",""),
               ifelse(show_progress,"--show-progress",""),
               "-P",output_path,
               ifelse(force_overwrite,"","--no-clobber"))
  # print(cmd)
  system(paste(cmd,"&& echo '+ wget download complete.'"))
  return(out.file)
}



#' axel
#'
#' R wrapper for axel, which enables multi-threaded download of a single large file.
#' @family download functions
#' @seealso \url{https://github.com/axel-download-accelerator/axel/}
#' @keywords internal
axel <- function(input_url,
                 output_path,
                 background=F,
                 nThreads=4,
                 force_overwrite=F,
                 quiet=F,
                 alternate=T){
  dir.create(output_path, showWarnings = F, recursive = T)
  out.file <- file.path(output_path,basename(input_url))
  if(force_overwrite){
    print("+ Overwriting pre-existing file.")
    suppressWarnings(file.remove(out.file))
  }

  cmd <- paste("axel",input_url,
               "-n",nThreads,
               ifelse(force_overwrite,"","--no-clobber"),
               "-o",out.file,
               ifelse(quiet,"-q",""),
               # ifelse(alternate,"-a",""),
               ifelse(background,"& bg","")
  )
  # print(cmd)
  system(cmd)
  return(out.file)
}
