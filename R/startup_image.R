#' Print the echolocatoR image to the R console.
#'
#' @family general
#' @examples
#' startup_image()
#' @keywords internal
startup_image  <- function(){
  library(dplyr)
  try({
    col.text <- function(txt){
      c(txt,"\n") %>%
        crayon::blurred() %>%
        crayon::bgBlack() %>%
        # crayon::col_align(align = "left") %>%
        crayon::cyan() %>%
        cat()
    }
    col.text("))))))))))>>))))))))))>  E c h o l o c a t o R  <((((((((((<<((((((((((")
    col.text("")
    col.text("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ V1.0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    col.text("~~~~~~~~~~~~~~~~~~~~~~ © 2019 - Brian M. Schilder ~~~~~~~~~~~~~~~~~~~~~")
    col.text("~Department of Neuroscience, Department of Genetics & Genomic Sciences~")
    col.text("~~~~~~~~~~Icahn School of Medicine at Mount Sinai, NY, NYC, USA~~~~~~~~")
    col.text("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    grid::grid.newpage()
    img <- png::readPNG("./echolocatoR/images/echo_logo.png")
    grid::grid.raster(img)
  })
}
