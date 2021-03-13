




#' Genome build liftover
#'
#' Transfer your genomic coordinates from one genome build to another.
#'
#' \code{xLiftOver} was extracted from the \href{http://xgr.r-forge.r-project.org}(XGR} package.
xLiftOver <- function (data.file,
                       format.file = c("data.frame", "bed", "chr:start-end",
                                     "GRanges"),
                       build.conversion = c(NA, "hg38.to.hg19", "hg19.to.hg38",
                                            "hg19.to.hg18", "hg18.to.hg38", "hg18.to.hg19"),
                       merged = T,
          verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
          guid = NULL)
{
  startT <- Sys.time()
  message(paste(c("Start at ", as.character(startT)), collapse = ""),
          appendLF = T)
  message("", appendLF = T)
  format.file <- match.arg(format.file)
  build.conversion <- match.arg(build.conversion)
  if (verbose) {
    now <- Sys.time()
    message(sprintf("First, import the files formatted as '%s' (%s) ...",
                    format.file, as.character(now)), appendLF = T)
  }
  if (verbose) {
    now <- Sys.time()
    message(sprintf("\timport the data file (%s) ...", as.character(now)),
            appendLF = T)
  }
  if (is.matrix(data.file) | is.data.frame(data.file) | class(data.file) ==
      "GRanges") {
    data <- data.file
  }
  else if (!is.null(data.file) & any(!is.na(data.file))) {
    if (length(data.file) == 1) {
      if (file.exists(data.file)) {
        data <- utils::read.delim(file = data.file, header = F,
                                  row.names = NULL, stringsAsFactors = F)
        data <- unique(data[, 1])
      }
      else {
        data <- data.file
      }
    }
    else {
      data <- data.file
    }
  }
  else {
    warning("The file 'data.file' must be provided!\n")
    return(NULL)
  }
  if (verbose) {
    now <- Sys.time()
    message(sprintf("Second, construct GenomicRanges object (%s) ...",
                    as.character(now)), appendLF = T)
  }
  if (format.file == "data.frame") {
    if (ncol(data) >= 3) {
      data <- data
    }
    else if (ncol(data) == 2) {
      data <- cbind(data, data[, 2])
    }
    else {
      stop("Your input 'data.file' is not as expected!\n")
    }
    ind <- suppressWarnings(which(!is.na(as.numeric(data[,
                                                         2])) & !is.na(as.numeric(data[, 3]))))
    data <- data[ind, ]
    dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[,
                                                                 1]), ranges = IRanges::IRanges(start = as.numeric(data[,
                                                                                                                        2]), end = as.numeric(data[, 3])), strand = S4Vectors::Rle(rep("*",
                                                                                                                                                                                       nrow(data))))
  }
  else if (format.file == "chr:start-end") {
    input <- do.call(rbind, strsplit(data[, 1], ":|-"))
    if (ncol(input) >= 3) {
      data <- input
    }
    else if (ncol(input) == 2) {
      data <- cbind(input, input[, 2])
    }
    else {
      stop("Your input 'data.file' does not meet the format 'chr:start-end'!\n")
    }
    ind <- suppressWarnings(which(!is.na(as.numeric(data[,
                                                         2])) & !is.na(as.numeric(data[, 3]))))
    data <- data[ind, ]
    dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[,
                                                                 1]), ranges = IRanges::IRanges(start = as.numeric(data[,
                                                                                                                        2]), end = as.numeric(data[, 3])), strand = S4Vectors::Rle(rep("*",
                                                                                                                                                                                       nrow(data))))
  }
  else if (format.file == "bed") {
    ind <- suppressWarnings(which(!is.na(as.numeric(data[,
                                                         2])) & !is.na(as.numeric(data[, 3]))))
    data <- data[ind, ]
    dGR <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data[,
                                                                 1]), ranges = IRanges::IRanges(start = as.numeric(data[,
                                                                                                                        2]) + 1, end = as.numeric(data[, 3])), strand = S4Vectors::Rle(rep("*",
                                                                                                                                                                                           nrow(data))))
  }
  else if (format.file == "GRanges") {
    dGR <- data
  }
  if (verbose) {
    now <- Sys.time()
    message(sprintf("Third, lift intervals between genome builds '%s' (%s) ...",
                    build.conversion, as.character(now)), appendLF = T)
  }
  chains <- xRDataLoader(RData.customised = "chain", RData.location = RData.location,
                         guid = guid, verbose = verbose)
  chain <- ""
  eval(parse(text = paste("chain <- chains$", build.conversion,
                          sep = "")))
  suppressMessages(res_GRL <- rtracklayer::liftOver(dGR, chain))
  res_GR <- BiocGenerics::unlist(res_GRL)
  if (merged) {
    mcols_data <- GenomicRanges::mcols(dGR)
    if (is.null(names(dGR))) {
      names(dGR) <- 1:length(dGR)
    }
    names_data <- names(dGR)
    if (verbose) {
      now <- Sys.time()
      message(sprintf("Finally, keep the first range if multiple found (%s) ...",
                      as.character(now)), appendLF = T)
    }
    res_df <- GenomicRanges::as.data.frame(res_GR, row.names = NULL)
    uid <- names(res_GR)
    res_ls <- split(x = res_df[, c(1:3, 5)], f = uid)
    ls_df <- lapply(res_ls, function(x) {
      c(as.character(unique(x$seqnames))[1], min(x$start),
        max(x$end), as.character(unique(x$strand))[1])
    })
    df <- do.call(rbind, ls_df)
    gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(df[,
                                                              1]), ranges = IRanges::IRanges(start = as.numeric(df[,
                                                                                                                   2]), end = as.numeric(df[, 3])), strand = S4Vectors::Rle(df[,
                                                                                                                                                                               4]))
    ind <- match(rownames(df), names_data)
    names(gr) <- names_data[ind]
    GenomicRanges::mcols(gr) <- mcols_data[ind, ]
    res_GR <- gr
  }
  endT <- Sys.time()
  message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
          appendLF = T)
  runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"),
                                 strptime(startT, "%Y-%m-%d %H:%M:%S"), units = "secs"))
  message(paste(c("Runtime in total is: ", runTime, " secs\n"),
                collapse = ""), appendLF = T)
  invisible(res_GR)
}
