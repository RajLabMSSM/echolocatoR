.onLoad <- function(libname, pkgname){
  .datatable.aware <- TRUE
}

.onAttach <- function(libname, pkgname){
  ## Suppress the decorative startup banner during R CMD check
  ## to avoid polluting check output with braille/unicode art.
  if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""))) return()
  msg_out <- utils::capture.output({
    msg_err <- utils::capture.output(
      startup(package = pkgname),
      type = "message"
    )
  })
  all_msg <- c(msg_out, msg_err)
  packageStartupMessage(paste(all_msg, collapse = "\n"))
}
