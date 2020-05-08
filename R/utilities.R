





#' printer
#'
#' Concatenate and print any number of items.
#'
#' @param ...
#' @return character
#' @examples
#' n.snps <- 50
#' printer("echolocatoR::","Processing",n.snps,"SNPs...")
printer <- function(..., v=T){if(v){print(paste(...))}}




#' Ensembl IDs -> Gene symbols
#'
#' Print the echolocatoR image to the R console.
#'
#' @examples
#' startup_image()
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




#' Bind stored files
#'
#' Rapidly read a list of files from storage and concatenate them by rows.
#'
#' @examples
#' file.list <- c("data/file1.tsv", "data/file2.tsv","new_data/file3/tsv")
#' dat <- rbind.file.list(file.list = file.list)
rbind.file.list <- function(file.list,
                            verbose=T,
                            nCores=4){
  merged.dat <- parallel::mclapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x)
    return(dat)
  }, mc.cores = nCores) %>% data.table::rbindlist(fill=T)
  return(merged.dat)
}




#' tryCatch extension
#'
#' Extension of tryCatch function.
#'
#' @family developer functions
#' @param input Function input.
#' @param func Function.
#' @examples
#' finemap_DT$effective_sample_size <- effective_sample_size(finemap_DT)
tryFunc <- function(input, func) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      func(input)
      # The return value of `readLines()` is the actual value
      # that will be returned in case there is no condition
      # (e.g. warning or error).
      # You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", input))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", input))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message(paste("Processed URL:", input))
      message("Some other message at the end")
    }
  )
  return(out)
}




#' Interactive DT
#'
#' Generate an interactive data table with download buttons.
#'
#' @param data.frame
#' @example
#' createDT(mtcars)
createDT <- function(DF, caption="", scrollY=400){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}




#' Interactive DT (html)
#'
#' Generate an interactive data table with download buttons.
#' Use this function when manually constructing rmarkdown chunks using cat() in a for loop.
#'
#' @param data.frame
#' @example
#' createDT_html(mtcars)
createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}




#' Replace items within DT object
#'
#' Annoyingly, there is no native function to do simple find-and-replace in the `DT` library.
#'
#' @param data.frame
#' @return data.frame
#' @example
#' dat <- dt.replace(mtcars, .1, NA)
dt.replace <- function(DT, target, replacement){
  for(col in names(DT)) set(DT, i=which(DT[[col]]==target), j=col, value=replacement)
  return(DT)
}


