#' Startup messages
#'
#' Startup messages with info on package name, version, citation, GitHub, etc.
#' @keywords internal
#' @import cli
startup <- function(package="echolocatoR",
                    add_batty=FALSE){

  ref <- gsub("\\[|\\*","",
              strsplit(citation("echolocatoR")$textVersion,"\\]")[[1]][1])
  indent <- 5
  exdent <- 5
  width <- (cli::console_width()*.75) - indent - exdent

  bat_icon <- function(n=1, sep="\ \ "){
    paste(rep("\U0001F987",n), collapse = sep)
  }

  #### Make waves function #####
  make_waves <- function(frequency=50,
                         width=cli::console_width(),
                         fun=sin,
                         cat_now=TRUE,
                         align="left"){
    obj <- cli::spark_line(scales::rescale(fun(seq(0, frequency, length = width))))
    obj <- cli::ansi_align(cli::col_br_cyan(obj),
                    align = align,
                    width = cli::console_width()*.95)
    if(cat_now){
      cli::cat_line(obj)
    } else {return(obj)}
  }

  #### Make padding ####
  add_pad <- function(txt,
                      pad=" ",
                  width=cli::console_width()#nchar(make_waves(cat_now = FALSE), type="width")
                  ){
    extra_space <- width - utf8::utf8_width(txt) #nchar(x, type="width")
    padded_txt <- stringr::str_pad(txt, width = extra_space,
                                   pad = pad, side = "right")
    return(padded_txt)
  }

  # border <- cli::col_br_cyan(cli::cli_h1(""))

  if(add_batty){
    batty <- readLines(
      system.file(package = "echolocatoR","extdata/ascii_batty.txt")
    )
    cat(cli::col_br_magenta(paste(batty, collapse = "\n")))
    cli::cat_line()
  }

  bat_title <- paste(bat_icon(3),
                     paste(strsplit(package,"")[[1]], collapse = " "),
                     bat_icon(3))
  {
      align <- "left"
      freqs <- c(50,100,200)
      #### Waves in ####
      for(x in freqs){
        make_waves(x, align=align)
      }
      #### Package name ####
      cli::cli_h1(cli::col_br_cyan(bat_title))
      #### Package version ####
      cli::cli_h1(cli::col_br_cyan(paste0("v",utils::packageVersion(package))))
    ##### Waves out ####
    for(x in rev(freqs)){
      make_waves(x,align=align)
    }
  }
  ##### Additional text #####
  {
    #### Citation
    cli::cat_line(
      cli::ansi_align(
        cli::col_cyan(paste0(cli::symbol$circle_circle,
                             " If you use ",package,", please cite:")
        ),align = "left"
      )
    )
    cli::cat_line(
      cli::ansi_align(
        cli::col_br_magenta(
          cli::ansi_strwrap(
            paste(cli::symbol$play,ref),
            width = width, indent = indent, exdent = exdent)
        ),align = "left"
      )
    )
    #### Issues ####
    cli::cat_line(
      cli::col_cyan(
        paste(cli::symbol$circle_circle,
              "Please report any bugs/feature requests on GitHub:")
      )
    )
    cli::cat_line(
      cli::col_br_magenta(
        cli::ansi_strwrap(
          paste(cli::symbol$play,
                cli::style_hyperlink(text = file.path("https://github.com/RajLabMSSM",package,"issues"),
                                     url = file.path("https://github.com/RajLabMSSM",package,"issues"))
                ) ,
          width = width, indent = indent, exdent = exdent)
      )
    )
    #### Contributions ####
    cli::cat_line(
      cli::col_cyan(
        paste(cli::symbol$circle_circle,"Contributions are welcome!:")
      )
    )
    cli::cat_line(
      cli::col_br_magenta(
        cli::ansi_strwrap(
          paste(cli::symbol$play,
                cli::style_hyperlink(text = file.path("https://github.com/RajLabMSSM",package,"pulls"),
                                     url = file.path("https://github.com/RajLabMSSM",package,"pulls"))

                ),
          width = width, indent = indent, exdent = exdent)
      )
    )
  }
  border <- cli::col_br_cyan(cli::cli_h1(""))

 # cli::ansi_with_hidden_cursor(cli::get_spinner(c("moon"))
 # cat(cli::ansi_html("<a href='https://doi.org/10.1093/bioinformatics/btab658'>link</a>)"))
}
