

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
                       alternate=T,
                       check_certificates=F,
                       verbose=T){
  if(download_method=="axel"){
    axel_avail <- length(system("which axel",intern = T))!=0
    if(axel_avail){
      out_file <- axel(input_url=input_url,
                       output_path=output_path,
                       background=background,
                       nThread=nThread,
                       force_overwrite=force_overwrite,
                       quiet=T, # output hella long otherwise...
                       alternate=alternate,
                       check_certificates=check_certificates)
    } else {
      printer("+ DOWNLOADER:: Axel not available. Defaulting to wget.",v=verbose);
      download_method <- "wget"
    }

  }
  if(download_method=="wget"){
    wget_avail <- length(system("which wget",intern = T))!=0
    if(wget_avail){
      out_file <- wget(input_url,
                       output_path,
                       background=background,
                       force_overwrite=force_overwrite,
                       quiet=quiet,
                       show_progress=show_progress,
                       continue=continue,
                       check_certificates=check_certificates)
    } else {
      stop("+ DOWNLOADER:: Please install wget or axel first.");
    }
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
                 continue=T,
                 check_certificates=F){
  # https://stackoverflow.com/questions/21365251/how-to-run-wget-in-background-for-an-unattended-download-of-files
  ## -bqc makes wget run in the background quietly
  dir.create(output_path, showWarnings = F, recursive = T)
  out_file <- file.path(output_path,basename(input_url))
  cmd <- paste("wget",input_url,
               "-np",
               ## Checking certificates can sometimes cause issues
               if(check_certificates) "" else "--no-check-certificate",
               if(background) "-b" else "",
               if(continue) "-c" else "",
               if(quiet) "-q" else "",
               if(show_progress) "--show-progress" else "",
               "-P",output_path,
               if(force_overwrite) "" else "--no-clobber"
  )
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
                 quiet=T,
                 alternate=T,
                 check_certificates=F){
  dir.create(output_path, showWarnings = F, recursive = T)
  out_file <- file.path(output_path,basename(input_url))
  if(force_overwrite){
    print("+ Overwriting pre-existing file.")
    suppressWarnings(file.remove(out_file))
  }
  cmd <- paste("axel",input_url,
               "-n",nThread,
               ## Checking certificates can sometimes cause issues
               if(check_certificates) "" else "--insecure",
               if(force_overwrite) "" else "--no-clobber",
               "-o",out_file,
               if(quiet) "-q" else "",
               # ifelse(alternate,"-a",""),
               if(background) "& bg" else ""
  )
  # print(cmd)
  system(cmd)
  return(out_file)
}
