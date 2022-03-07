.onLoad <- function(libname, pkgname){
  banner <- paste(rep("~",100), collapse = "")
  indent <- ">>>>    "
  ref <- gsub("\\[|\\*","",
              strsplit(citation("echolocatoR")$textVersion,"\\]")[[1]][1])
  batty <- readLines(
    system.file(package = "echolocatoR","extdata/ascii_batty.txt")
  )

  txt <- paste(
    banner,
    cat(batty, sep="\n"),
    paste("🦇🦇🦇 )))))))))))>",
                "e c h o l o c a t o R",
                paste0("(v",utils::packageVersion("echolocatoR"),")"),
          "<((((((((((( 🦇🦇🦇"),
    "",
    "If you use echolocatoR, please cite:",
    paste(indent,paste(strwrap(ref, 100), collapse = "\n")),
    "",
    "Please report any bugs/requests on GitHub:",
    paste(indent,"https://github.com/RajLabMSSM/echolocatoR/issues"),
    "",
    "Contributions are welcome!:",
    paste(indent,"https://github.com/RajLabMSSM/echolocatoR/pulls"),
    banner,
    sep = "\n"
   )
  message(txt)
}
