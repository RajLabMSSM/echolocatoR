
# R wrapper for faster downloading software


wget <- function(input.url,
                 output.path,
                 background=T,
                 force_overwrite=F,
                 quiet=F,
                 show_progress=T,
                 continue=T){
  # https://stackoverflow.com/questions/21365251/how-to-run-wget-in-background-for-an-unattended-download-of-files
  ## -bqc makes wget run in the background quietly
  dir.create(output.path, showWarnings = F, recursive = T)
  out.file <- file.path(output.path,basename(input.url))
  cmd <- paste("wget",input.url,
               "-np",
               ifelse(background,"-b",""),
               ifelse(continue,"-c",""),
               ifelse(quiet,"-q",""),
               ifelse(show_progress,"--show-progress",""),
               "-P",output.path,
               ifelse(force_overwrite,"","--no-clobber"))
  # print(cmd)
  system(paste(cmd,"&& echo '+ wget download complete.'"))
  return(out.file)
}

axel <- function(input.url,
                 output.path,
                 background=F,
                 nThreads=4,
                 force_overwrite=F,
                 quiet=F,
                 alternate=T){
  # https://github.com/axel-download-accelerator/axel/
  dir.create(output.path, showWarnings = F, recursive = T)
  out.file <- file.path(output.path,basename(input.url))
  if(force_overwrite){
    print("+ Overwriting pre-existing file.")
    suppressWarnings(file.remove(out.file))
  }

  cmd <- paste("axel",input.url,
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
