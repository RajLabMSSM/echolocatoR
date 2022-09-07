#' Startup messages
#'
#' Startup messages with info on package name, version, citation, GitHub, etc.
#' @param package R package name.
#' @param add_batty Add an ASCII bat.
#' @param color1 First color to use in the palette.
#' @param color2 Second color to use in the palette.
#' @param color3 Second color to use in the palette.
#' @returns Null
#'
#' @keywords internal
#' @import cli
#' @importFrom utils citation
#' @importFrom scales rescale
#' @importFrom utf8 utf8_width
startup <- function(package="echolocatoR",
                    add_batty=FALSE,
                    color1=cli::col_br_cyan,
                    color2=cli::col_br_magenta,
                    color3=cli::col_cyan){

  ref <- gsub("\\[|\\*","",
              strsplit(
                utils::citation("echolocatoR")$textVersion,"\\]"
                )[[1]][1])
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
    obj <- cli::spark_line(
      scales::rescale(fun(seq(0, frequency, length = width)))
      )
    obj <- cli::ansi_align(color1(obj),
                    align = align,
                    width = cli::console_width()*.95)
    if(cat_now){
      cli::cat_line(obj)
    } else {return(obj)}
  }

  #### Make padding ####
  add_pad <- function(txt,
                      pad=" ",
                      #nchar(make_waves(cat_now = FALSE), type="width")
                  width=cli::console_width()
                  ){
    extra_space <- width - utf8::utf8_width(txt) #nchar(x, type="width")
    padded_txt <- stringr::str_pad(txt, width = extra_space,
                                   pad = pad, side = "right")
    return(padded_txt)
  }

  # border <- color1(cli::cli_h1(""))

  if(add_batty){
    batty <- readLines(
      system.file(package = "echolocatoR","extdata/ascii_batty.txt")
    )
    cat(color2(paste(batty, collapse = "\n")))
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
      cli::cli_h1(color1(bat_title))
      #### Package version ####
      cli::cli_h1(color1(paste0("v",utils::packageVersion(package))))
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
        color3(paste0(cli::symbol$circle_circle,
                      " If you use ",cli::style_bold(package),
                      " or any of the ",cli::style_italic("echoverse"),
                      " subpackages, please cite:")
        ),align = "left"
      )
    )
    cli::cat_line(
      cli::ansi_align(
        color2(
          cli::ansi_strwrap(
            paste(cli::symbol$play,ref),
            width = width, indent = indent, exdent = exdent)
        ),align = "left"
      )
    )
    #### Issues ####
    cli::cat_line(
      color3(
        paste(cli::symbol$circle_circle,
              "Please report any bugs/feature requests on GitHub:")
      )
    )
    cli::cat_line(
      color2(
        cli::ansi_strwrap(
          paste(cli::symbol$play,
                cli::style_hyperlink(
                  text = file.path(
                    "https://github.com/RajLabMSSM",package,"issues"
                    ),
                  url = file.path(
                    "https://github.com/RajLabMSSM",package,"issues")
                  )
                ) ,
          width = width, indent = indent, exdent = exdent)
      )
    )
    #### Contributions ####
    cli::cat_line(
      color3(
        paste(cli::symbol$circle_circle,"Contributions are welcome!:")
      )
    )
    cli::cat_line(
      color2(
        cli::ansi_strwrap(
          paste(cli::symbol$play,
                cli::style_hyperlink(
                  text = file.path(
                    "https://github.com/RajLabMSSM",package,"pulls"
                    ),
                  url = file.path(
                    "https://github.com/RajLabMSSM",package,"pulls")
                  )
                ),
          width = width, indent = indent, exdent = exdent)
      )
    )
  }
  border <- color1(cli::cli_h1(""))

 # cli::ansi_with_hidden_cursor(cli::get_spinner(c("moon"))
 # cat(cli::ansi_html("<a href='https://doi.org/10.1093/bioinformatics/btab658'>link</a>)"))
}
