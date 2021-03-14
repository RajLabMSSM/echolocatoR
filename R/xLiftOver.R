




#' Genome build liftover
#'
#' Transfer your genomic coordinates from one genome build to another.
#'
#' \code{xLiftOver} was extracted from the \href{http://xgr.r-forge.r-project.org}(XGR} package.
#' @export
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







xRDataLoader <- function (RData = c(NA, "GWAS2EF", "GWAS_LD", "IlluminaHumanHT",
                    "IlluminaOmniExpress", "ig.DO", "ig.EF", "ig.GOBP", "ig.GOCC",
                    "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP",
                    "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP",
                    "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA",
                    "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1",
                    "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall",
                    "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME",
                    "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN",
                    "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC",
                    "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7",
                    "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.egPfam",
                    "org.Hs.string", "org.Hs.PCommons_DN", "org.Hs.PCommons_UN"),
          RData.customised = NULL, verbose = T, RData.location = "http://galahad.well.ox.ac.uk/bigdata",
          guid = NULL)
{
  startT <- Sys.time()
  if (verbose) {
    message(paste(c("Start at ", as.character(startT)), collapse = ""),
            appendLF = TRUE)
    message("", appendLF = TRUE)
  }
  RData <- RData[1]
  if (is.na(RData) & !is.null(RData.customised)) {
    RData <- RData.customised
  }
  else if (is.na(RData) & is.null(RData.customised)) {
    stop("There is no input! Please input one of two parameters ('RData' or 'RData.customised').\n")
  }
  RData <- gsub(".RData$", "", RData, ignore.case = T, perl = T)
  RData <- gsub(".RDa$", "", RData, ignore.case = T, perl = T)
  flag_osf <- F
  if (!is.null(guid) && nchar(guid) == 5) {
    pkgs <- c("osfr")
    if (all(pkgs %in% rownames(utils::installed.packages()))) {
      tmp <- sapply(pkgs, function(pkg) {
        requireNamespace(pkg, quietly = T)
      })
      if (all(tmp)) {
        prj <- fls <- res <- NULL
        if (all(class(suppressWarnings(try(eval(parse(text = paste0("prj<-osfr::osf_retrieve_node(guid)"))),
                                           T))) != "try-error")) {
          target <- paste0(RData, ".RData")
          eval(parse(text = noquote(paste0("fls <- osfr::osf_ls_files(prj, type=\"file\", pattern=target, n_max=Inf)"))))
          if (nrow(fls) > 0) {
            ind <- match(fls$name, target)
            ind <- ind[!is.na(ind)]
            if (length(ind) == 1) {
              fl <- fls[ind, ]
              destfile <- file.path(tempdir(), fl$name)
              eval(parse(text = paste0("res <- fl %>% osfr::osf_download(overwrite=T, path=destfile)")))
              if (file.exists(res$local_path)) {
                out <- get(load(res$local_path))
                load_RData <- sprintf("'%s' at %s", prj$name,
                                      paste0("https://osf.io/", prj$id))
                RData <- target
                flag_osf <- T
              }
            }
          }
        }
      }
    }
  }
  my_https_downloader <- function(url, method = c("auto", "internal",
                                                  "wininet", "libcurl", "wget", "curl"), quiet = T, mode = c("w",
                                                                                                             "wb", "a", "ab"), cacheOK = T, extra = getOption("download.file.extra")) {
    method <- match.arg(method)
    mode <- match.arg(mode)
    tdir <- tempdir()
    destfile <- file.path(tdir, "temp.RData")
    unlink(destfile, recursive = T, force = T)
    if (base::grepl("^https?://", url)) {
      isR32 <- base::getRversion() >= "3.2"
      if (.Platform$OS.type == "windows") {
        if (isR32) {
          method <- "wininet"
        }
        else {
          seti2 <- utils::"setInternet2"
          internet2_start <- seti2(NA)
          if (!internet2_start) {
            on.exit(suppressWarnings(seti2(internet2_start)))
            suppressWarnings(seti2(TRUE))
          }
          method <- "internal"
        }
      }
      else {
        if (isR32 && capabilities("libcurl")) {
          method <- "libcurl"
        }
        else if (nzchar(Sys.which("wget")[1])) {
          method <- "wget"
        }
        else if (nzchar(Sys.which("curl")[1])) {
          method <- "curl"
          orig_extra_options <- getOption("download.file.extra")
          on.exit(options(download.file.extra = orig_extra_options))
          options(download.file.extra = paste("-L", orig_extra_options))
        }
        else if (nzchar(Sys.which("lynx")[1])) {
          method <- "lynx"
        }
        else {
          stop("no download method found")
        }
      }
    }
    else {
    }
    if (class(suppressWarnings(try(utils::download.file(url,
                                                        destfile = destfile, method = method, quiet = quiet,
                                                        mode = mode, cacheOK = cacheOK, extra = extra), T))) ==
        "try-error") {
      res_RData <- NULL
      res_flag <- F
    }
    if (file.exists(destfile) & file.info(destfile)$size !=
        0) {
      if (class(suppressWarnings(try(load(destfile), T))) ==
          "try-error") {
        res_RData <- NULL
        res_flag <- F
      }
      else {
        res_RData <- get(load(destfile))
        res_flag <- T
      }
    }
    else {
      res_RData <- NULL
      res_flag <- F
    }
    res <- list(RData = res_RData, flag = res_flag)
    invisible(res)
  }
  if (!flag_osf) {
    path_host <- gsub("/$", "", RData.location)
    if (path_host == "" || length(path_host) == 0 || is.na(path_host)) {
      path_host <- "https://github.com/hfang-bristol/RDataCentre/blob/master/Portal"
    }
    load_remote <- paste(path_host, "/", RData, ".RData",
                         sep = "")
    load_local1 <- file.path(path_host, paste("data/", RData,
                                              ".RData", sep = ""))
    load_local2 <- file.path(path_host, paste(RData, ".RData",
                                              sep = ""))
    load_package <- RData
    if (1) {
      RData_local <- c(load_local1, load_local2)
      load_flag <- sapply(RData_local, function(x) {
        if (.Platform$OS.type == "windows")
          x <- gsub("/", "\\\\", x)
        ifelse(file.exists(x), TRUE, FALSE)
      })
      if (sum(load_flag) == 0) {
        flag_failed <- F
        if (length(grep("^https", load_remote, perl = T))) {
          if (length(grep("github", load_remote, perl = T))) {
            load_remote <- paste(load_remote, "?raw=true",
                                 sep = "")
          }
          res <- my_https_downloader(load_remote, mode = "wb")
          if (res$flag == F) {
            flag_failed <- T
          }
          else {
            eval(parse(text = paste(RData, " <- res$RData",
                                    sep = "")))
          }
        }
        else {
          res <- my_https_downloader(load_remote, mode = "wb")
          if (res$flag == F) {
            flag_failed <- T
          }
          else {
            eval(parse(text = paste(RData, " <- res$RData",
                                    sep = "")))
          }
        }
        if (flag_failed) {
          load_remotes <- c(paste("https://github.com/hfang-bristol/RDataCentre/blob/master/Portal/",
                                  RData, ".RData?raw=true", sep = ""), paste("http://galahad.well.ox.ac.uk/bigdata/",
                                                                             RData, ".RData", sep = ""), paste("http://galahad.well.ox.ac.uk/bigdata/",
                                                                                                               RData, ".RData", sep = ""))
          for (i in 1:length(load_remotes)) {
            load_remote <- load_remotes[i]
            if (verbose) {
              now <- Sys.time()
              message(sprintf("Attempt to download from %s (at %s)",
                              load_remote, as.character(now)), appendLF = T)
            }
            res <- my_https_downloader(load_remote, mode = "wb")
            if (res$flag == T) {
              break
            }
          }
          if (res$flag == F) {
            warnings("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")
            eval(parse(text = paste(RData, " <- res$RData",
                                    sep = "")))
          }
          else {
            eval(parse(text = paste(RData, " <- res$RData",
                                    sep = "")))
          }
        }
        load_RData <- load_remote
        out <- base::get(RData)
      }
      else {
        load_RData <- RData_local[load_flag]
        out <- base::get(load(load_RData))
      }
    }
    else {
      load_RData <- sprintf("package 'XGR' version %s",
                            utils::packageVersion("XGR"))
      out <- base::get(RData)
    }
  }
  if (verbose) {
    now <- Sys.time()
    if (!is.null(out)) {
      message(sprintf("'%s' (from %s) has been loaded into the working environment (at %s)",
                      RData, load_RData, as.character(now)), appendLF = T)
    }
    else {
      message(sprintf("'%s' CANNOT be loaded (at %s)",
                      RData, as.character(now)), appendLF = T)
    }
  }
  endT <- Sys.time()
  runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"),
                                 strptime(startT, "%Y-%m-%d %H:%M:%S"), units = "secs"))
  if (verbose) {
    message(paste(c("\nEnd at ", as.character(endT)), collapse = ""),
            appendLF = TRUE)
    message(paste(c("Runtime in total is: ", runTime, " secs\n"),
                  collapse = ""), appendLF = TRUE)
  }
  invisible(out)
}


