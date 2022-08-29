steps <- function(name = c("query",
                           "ld",
                           "filter",
                           "finemap",
                           "plot",
                           "postprocess")){
  name <- tolower(name)
  div <- paste(rep(cli::symbol$play,3),collapse = "")
  if("query" %in% name){
    cli::cli_h1("")
    cli::cli_h1(paste("Step 1",div,"Query \U0001f50e"))
    cli::cli_h1("")
  }
  if("ld" %in% name){
    cli::cli_h1("")
    cli::cli_h1(paste("Step 2",div,"Extract Linkage Disequilibrium \U0001f517"))
    cli::cli_h1("")
  }
  if("filter" %in% name){
    cli::cli_h1("")
    cli::cli_h1(paste("Step 3",div,"Filter SNPs \U0001f6b0"))
    cli::cli_h1("")
  }
  if("finemap" %in% name){
    cli::cli_h1("")
    cli::cli_h1(paste("Step 4",div,"Fine-map \U0001f50a"))
    cli::cli_h1("")
  }
  if("plot" %in% name){
    cli::cli_h1("")
    cli::cli_h1(paste("Step 5",div,"Plot \U0001f4c8"))
    cli::cli_h1("")
  }
  if("postprocess" %in% name){
    cli::cli_h1("")
    cli::cli_h1(paste("Step 6",div,"Postprocess data \U0001f381"))
    cli::cli_h1("")
  }
}
