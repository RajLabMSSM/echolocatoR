report_time <- function(t1,
                        t2=NULL,
                        prefix="Done in:"){
  if(is.null(t2)) t2 <- Sys.time()
  out <- as.numeric(round(difftime(t2,t1,units = "m"),2))
  cat(
    cli::col_br_cyan(paste(
      prefix,cli::col_br_white(out," min")
    ),"\n"
    )
  )
}
