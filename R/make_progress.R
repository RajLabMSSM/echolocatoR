
make_progress <- function(X,
                          FUN=function(l)Sys.sleep(5/100),
                          ...,
                          apply_func=lapply,
                          total = length(X),
                          name=NULL,
                          .envir=parent.frame(),
                          cli.progress_show_after=0,
                          style=c("bar","squares","dot","fillsquares","classic"),
                          type=c("iterator", "tasks", "download", "custom"),
                          clear=FALSE) {
  #### Set up styles ####
  # ?cli::cli_progress_styles()
  options(cli.progress_show_after=cli.progress_show_after,
    cli.progress_bar_style = list(
    complete = cli::col_br_cyan(""),
    incomplete = cli::col_br_cyan("))>")
    ),
    cli.spinner = "dots" # "arrow3"
  )
  # X <- 1:100
  # names(X) <- make.unique(rep("a",100))
  # FUN=function(l)Sys.sleep(5/100)
  f <- function(X, FUN, name, clear=TRUE, ...){
    lapply(cli::cli_progress_along(
      x = X,
      name = name,
      format = paste0("{cli::col_br_white(cli::pb_spin)}",
                      " {cli::bg_white(cli::col_br_magenta(cli::pb_name))}",
                      " {names(X)[cli::pb_current]}",
                      " {paste0(cli::col_cyan('('),cli::col_cyan(cli::pb_percent),cli::col_cyan(')'))}",
                   " \U0001F9EC",
                   " {cli::pb_bar}",
                   "\U0001F987"),
      clear = clear),
      FUN,
      ...)
  }
  res <- f(X = X, FUN = FUN, name = name, clear = clear)
  #### There's also pbmcapply
  # res <- pbmcapply::pbmclapply(1:10, function(x){Sys.sleep(5/100)}, )
  return(res)
}
