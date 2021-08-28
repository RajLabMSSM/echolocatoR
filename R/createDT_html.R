#' Interactive DT (html)
#'
#' Generate an interactive data table with download buttons.
#' Use this function when manually constructing rmarkdown chunks using cat() in a for loop.
#'
#' @family general
#' @keywords internal
createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}
