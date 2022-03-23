
make_progress <- function(X,
                          FUN,
                          ...,
                          apply_func=lapply,
                          total = length(X),
                          name="Processing",
                          .envir=parent.frame(),
                          cli.progress_show_after=0,
                          style=c("bar","squares","dot","fillsquares","bar"),
                          type=c("iterator", "tasks", "download", "custom"),
                          clear=FALSE) {
  #### Set up styles ####
  cli::cli_progress_styles()
  options(cli.progress_show_after = cli.progress_show_after,
          cli.progress_bar_style = style[1]
  )
  cli::cli_progress_bar(name = name,
                        total = total,
                        type = type[1],
                        clear = clear)
  # for (i in seq_len(total)) {
  #   Sys.sleep(5/100)
  #   cli::cli_progress_update()
  # }
  # make_progress()
  #### Experimental version #####
  res <- apply_func(cli::cli_progress_along(x = X,
                                            name = name,
                                            total = total,
                                            .envir = .envir),
                    FUN, ...)
  #### There's also pbmcapply
  # res <- pbmcapply::pbmclapply(1:10, function(x){Sys.sleep(5/100)}, )
  return(res)
}
