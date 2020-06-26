

#' Downloaders
#'
#' R wrapper for wget and axel
#' @family downloaders
#' @keywords internal
downloader <- function(input_url,
                       output_path,
                       download_method="axel",
                       background=F,
                       force_overwrite=F,
                       quiet=F,
                       show_progress=T,
                       continue=T,

                       nThread=4,
                       alternate=T){
  if(download_method=="axel"){
    out_file <- axel(input_url=input_url,
                     output_path=output_path,
                     background=background,
                     nThread=nThread,
                     force_overwrite=force_overwrite,
                     quiet=quiet,
                     alternate=alternate)
  }
  if(download_method=="wget"){
    out_file <- wget(input_url,
                     output_path,
                     background=background,
                     force_overwrite=force_overwrite,
                     quiet=quiet,
                     show_progress=show_progress,
                     continue=continue)
  }
  return(out_file)
}




#' wget
#'
#' R wrapper for wget
#' @family downloaders
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
  out_file <- file.path(output_path,basename(input_url))
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
  return(out_file)
}



#' axel
#'
#' R wrapper for axel, which enables multi-threaded download of a single large file.
#' @family downloaders
#' @seealso \url{https://github.com/axel-download-accelerator/axel/}
#' @keywords internal
axel <- function(input_url,
                 output_path,
                 background=F,
                 nThread=4,
                 force_overwrite=F,
                 quiet=F,
                 alternate=T){
  dir.create(output_path, showWarnings = F, recursive = T)
  out_file <- file.path(output_path,basename(input_url))
  if(force_overwrite){
    print("+ Overwriting pre-existing file.")
    suppressWarnings(file.remove(out_file))
  }
  cmd <- paste("axel",input_url,
               "-n",nThread,
               ifelse(force_overwrite,"","--no-clobber"),
               "-o",out_file,
               ifelse(quiet,"-q",""),
               # ifelse(alternate,"-a",""),
               ifelse(background,"& bg","")
  )
  # print(cmd)
  system(cmd)
  return(out_file)
}
