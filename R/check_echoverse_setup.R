#' Check echoverse setup
#'
#' Diagnose the installation status of all echoverse packages,
#' system dependencies, external tools, and Python/conda environments.
#' Prints an actionable report with platform-specific fix commands.
#'
#' @param echoverse_pkgs Character vector of echoverse package names to check.
#' @param key_deps Character vector of key non-echoverse dependency
#'   package names to check.
#' @param verbose Print detailed results. Default \code{TRUE}.
#' @returns A list (invisibly) with elements:
#' \describe{
#'   \item{packages}{data.frame of echoverse package status}
#'   \item{system}{list of system dependency checks}
#'   \item{tools}{list of external tool checks}
#'   \item{python}{list of Python/conda checks}
#'   \item{pass}{logical, TRUE if all critical checks pass}
#' }
#' @export
#' @examples
#' \dontrun{
#' ## Full diagnostic report
#' results <- check_echoverse_setup()
#'
#' ## Quick check with a subset of packages
#' results <- check_echoverse_setup(
#'     echoverse_pkgs = c("echodata"),
#'     key_deps = c("data.table"),
#'     verbose = FALSE
#' )
#' }
check_echoverse_setup <- function(
    echoverse_pkgs = c(
        "echolocatoR", "echodata", "echotabix", "echoannot",
        "echoconda", "echoLD", "echoplot", "echofinemap",
        "catalogueR", "downloadR", "echogithub", "devoptera",
        "echodeps", "echoAI", "echoverseTemplate"
    ),
    key_deps = c(
        "data.table", "BiocManager", "reticulate",
        "susieR", "coloc", "MungeSumstats",
        "VariantAnnotation", "rtracklayer", "GenomicRanges",
        "ggbio", "basilisk", "Rsamtools"
    ),
    verbose = TRUE) {

    os <- get_os()
    results <- list(pass = TRUE)

    ## ---- Header ----
    if (verbose) {
        cli::cli_h1("echoverse setup check")
        cli::cli_text("OS: {os} | R: {R.version.string}")
        cat("\n")
    }

    ## ---- 1. echoverse packages ----
    pkg_status <- check_packages(echoverse_pkgs, verbose = verbose)
    results$packages <- pkg_status

    ## ---- 2. Key non-echoverse dependencies ----
    dep_status <- check_packages(key_deps, label = "Key dependencies",
                                 verbose = verbose)
    results$dependencies <- dep_status

    ## ---- 3. susieR version ----
    if (verbose) cli::cli_h2("Version checks")
    susie_ok <- check_pkg_version("susieR", "0.12.0", verbose = verbose)
    if (!susie_ok) results$pass <- FALSE

    ## ---- 4. System dependencies ----
    if (verbose) cli::cli_h2("System dependencies")
    sys_results <- check_system_deps(os = os, verbose = verbose)
    results$system <- sys_results
    if (!sys_results$all_ok) results$pass <- FALSE

    ## ---- 5. Python / conda ----
    if (verbose) cli::cli_h2("Python / conda")
    py_results <- check_python_setup(verbose = verbose)
    results$python <- py_results

    ## ---- 6. External tools ----
    if (verbose) cli::cli_h2("External fine-mapping tools")
    tool_results <- check_external_tools(verbose = verbose)
    results$tools <- tool_results

    ## ---- 7. GITHUB_TOKEN ----
    if (verbose) cli::cli_h2("GitHub authentication")
    gh_token <- check_github_token(verbose = verbose)
    results$github_token <- gh_token

    ## ---- Summary ----
    if (verbose) {
        cat("\n")
        cli::cli_h1("Summary")
        n_pkg <- sum(pkg_status$loadable)
        n_total <- nrow(pkg_status)
        n_dep <- sum(dep_status$loadable)
        n_dep_total <- nrow(dep_status)
        cli::cli_text(
            "echoverse packages: {n_pkg}/{n_total} loadable"
        )
        cli::cli_text(
            "Key dependencies: {n_dep}/{n_dep_total} loadable"
        )
        if (results$pass) {
            cli::cli_alert_success("All critical checks passed.")
        } else {
            cli::cli_alert_danger(
                "Some checks failed. See messages above for fixes."
            )
        }
    }

    invisible(results)
}


# ---- Internal helpers ----

#' @keywords internal
get_os <- function() {
    sys <- Sys.info()["sysname"]
    switch(tolower(sys),
           "darwin" = "macOS",
           "linux" = "Linux",
           "windows" = "Windows",
           as.character(sys))
}

#' @keywords internal
check_packages <- function(pkgs, label = NULL, verbose = TRUE) {
    if (is.null(label)) label <- "echoverse packages"
    if (verbose) cli::cli_h2(label)

    status <- data.frame(
        package = pkgs,
        installed = logical(length(pkgs)),
        loadable = logical(length(pkgs)),
        version = character(length(pkgs)),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(pkgs)) {
        pkg <- pkgs[i]
        status$installed[i] <- requireNamespace(pkg, quietly = TRUE)
        if (status$installed[i]) {
            status$version[i] <- tryCatch(
                as.character(utils::packageVersion(pkg)),
                error = function(e) "?"
            )
            status$loadable[i] <- tryCatch({
                suppressPackageStartupMessages(
                    loadNamespace(pkg)
                )
                TRUE
            }, error = function(e) FALSE)
        }
        if (verbose) {
            if (status$loadable[i]) {
                cli::cli_alert_success("{pkg} ({status$version[i]})")
            } else if (status$installed[i]) {
                cli::cli_alert_warning(
                    "{pkg} ({status$version[i]}) - installed but fails to load"
                )
            } else {
                cli::cli_alert_danger("{pkg} - not installed")
            }
        }
    }
    status
}

#' @keywords internal
check_pkg_version <- function(pkg, min_version, verbose = TRUE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (verbose) {
            cli::cli_alert_danger(
                "{pkg} not installed (need >= {min_version})"
            )
            cli::cli_text(
                '
    Fix: {.code install.packages("{pkg}")}'
            )
        }
        return(FALSE)
    }
    ver <- utils::packageVersion(pkg)
    ok <- ver >= min_version
    if (verbose) {
        if (ok) {
            cli::cli_alert_success("{pkg} {ver} (>= {min_version})")
        } else {
            cli::cli_alert_danger(
                "{pkg} {ver} is too old (need >= {min_version})"
            )
            cli::cli_text(
                '
    Fix: {.code devtools::install_github("stephenslab/susieR")}'
            )
        }
    }
    ok
}

#' @keywords internal
check_system_deps <- function(os, verbose = TRUE) {
    results <- list(all_ok = TRUE)

    ## libxml2
    xml_ok <- requireNamespace("XML", quietly = TRUE) ||
              requireNamespace("xml2", quietly = TRUE)
    results$libxml2 <- xml_ok
    if (verbose) {
        if (xml_ok) {
            cli::cli_alert_success("XML/xml2 package available (libxml2 OK)")
        } else {
            cli::cli_alert_danger("XML/xml2 not available (libxml2 missing)")
            results$all_ok <- FALSE
            if (os == "macOS") {
                cli::cli_text("
    Fix: {.code brew install libxml2}")
            } else if (os == "Linux") {
                cli::cli_text(
                    "
    Fix: {.code sudo apt-get install libxml2-dev}"
                )
            }
        }
    }

    ## zlib (needed by data.table fwrite compression)
    dt_ok <- tryCatch({
        tmp <- tempfile()
        data.table::fwrite(data.table::data.table(x = 1), tmp)
        unlink(tmp)
        TRUE
    }, error = function(e) FALSE)
    results$zlib <- dt_ok
    if (verbose) {
        if (dt_ok) {
            cli::cli_alert_success("data.table compression works (zlib OK)")
        } else {
            cli::cli_alert_warning(
                "data.table fwrite compression may be broken (zlib headers missing at compile time)"
            )
        }
    }

    ## Fortran compiler (needed by DescTools and others)
    fortran_ok <- requireNamespace("DescTools", quietly = TRUE)
    results$fortran <- fortran_ok
    if (verbose) {
        if (fortran_ok) {
            cli::cli_alert_success(
                "DescTools loadable (Fortran compiler was available)"
            )
        } else {
            cli::cli_alert_warning(
                "DescTools not available (may need Fortran compiler)"
            )
            if (os == "macOS") {
                cli::cli_text(
                    "
    Fix: {.code brew install gcc} (provides gfortran)"
                )
            } else if (os == "Linux") {
                cli::cli_text(
                    "
    Fix: {.code sudo apt-get install gfortran}"
                )
            }
        }
    }

    results
}

#' @keywords internal
check_python_setup <- function(verbose = TRUE) {
    results <- list()

    ## Python
    py <- Sys.which("python3")
    if (nchar(py) == 0) py <- Sys.which("python")
    results$python_found <- nchar(py) > 0
    if (verbose) {
        if (results$python_found) {
            py_ver <- tryCatch(
                system2(py, "--version", stdout = TRUE, stderr = TRUE),
                error = function(e) "unknown"
            )
            cli::cli_alert_success("Python: {py_ver} ({py})")
        } else {
            cli::cli_alert_danger("Python not found on PATH")
            cli::cli_text(
                "
    Fix: Install Python >= 3.7 from https://www.python.org"
            )
        }
    }

    ## conda
    conda <- ""
    if (requireNamespace("reticulate", quietly = TRUE)) {
        conda <- tryCatch(
            reticulate::conda_binary(),
            error = function(e) ""
        )
    }
    if (nchar(conda) == 0) conda <- Sys.which("conda")
    results$conda_found <- nchar(conda) > 0 && file.exists(conda)
    if (verbose) {
        if (results$conda_found) {
            cli::cli_alert_success("conda: {conda}")
        } else {
            cli::cli_alert_warning(
                "conda not found (needed for PolyFun, some LD methods)"
            )
            cli::cli_text(
                "
    Fix: Install Miniforge: https://github.com/conda-forge/miniforge"
            )
        }
    }

    ## echoR_mini env
    if (results$conda_found &&
        requireNamespace("reticulate", quietly = TRUE)) {
        envs <- tryCatch(
            reticulate::conda_list(),
            error = function(e) data.frame(name = character(0))
        )
        results$echoR_mini <- "echoR_mini" %in% envs$name
        if (verbose) {
            if (results$echoR_mini) {
                cli::cli_alert_success("conda env 'echoR_mini' exists")
            } else {
                cli::cli_alert_warning(
                    "conda env 'echoR_mini' not found (will be created on first use)"
                )
            }
        }
    }

    results
}

#' @keywords internal
check_external_tools <- function(verbose = TRUE) {
    results <- list()

    ## FINEMAP
    finemap_found <- tryCatch({
        exe <- echofinemap:::FINEMAP_find_executable(verbose = FALSE)
        file.exists(exe)
    }, error = function(e) FALSE)
    results$FINEMAP <- finemap_found
    if (verbose) {
        if (finemap_found) {
            cli::cli_alert_success("FINEMAP executable found")
        } else {
            cli::cli_alert_info(
                "FINEMAP not found (will be auto-downloaded on first use)"
            )
        }
    }

    ## PAINTOR
    paintor_found <- tryCatch({
        exe <- system.file("tools", "PAINTOR_V3.0", "PAINTOR",
                           package = "echofinemap")
        nchar(exe) > 0 && file.exists(exe)
    }, error = function(e) FALSE)
    results$PAINTOR <- paintor_found
    if (verbose) {
        if (paintor_found) {
            cli::cli_alert_success("PAINTOR executable found")
        } else {
            cli::cli_alert_info(
                "PAINTOR not compiled (requires manual build from echofinemap submodule)"
            )
        }
    }

    ## POLYFUN
    polyfun_found <- tryCatch({
        pf <- system.file("tools", "polyfun", package = "echofinemap")
        nchar(pf) > 0 && file.exists(file.path(pf, "polyfun.py"))
    }, error = function(e) FALSE)
    results$POLYFUN <- polyfun_found
    if (verbose) {
        if (polyfun_found) {
            cli::cli_alert_success("PolyFun found")
        } else {
            cli::cli_alert_info(
                "PolyFun not found (git submodule not checked out)"
            )
        }
    }

    ## genetics.binaRies (plink, gcta, bcftools)
    gbin_found <- requireNamespace("genetics.binaRies", quietly = TRUE)
    results$genetics_binaries <- gbin_found
    if (verbose) {
        if (gbin_found) {
            cli::cli_alert_success("genetics.binaRies installed (plink, gcta)")
        } else {
            cli::cli_alert_info(
                "genetics.binaRies not installed (needed for COJO methods)"
            )
            cli::cli_text(
                '
    Fix: {.code remotes::install_github("MRCIEU/genetics.binaRies")}'
            )
        }
    }

    if (verbose) {
        cli::cli_text("")
        cli::cli_alert_info(
            paste0(
                "FINEMAP, PAINTOR, and PolyFun are only needed if you ",
                "select those methods in finemap_loci()."
            )
        )
    }

    results
}

#' @keywords internal
check_github_token <- function(verbose = TRUE) {
    token <- Sys.getenv("GITHUB_TOKEN",
                 Sys.getenv("GITHUB_PAT",
                 Sys.getenv("GH_TOKEN", "")))
    ## Also check git credential store
    if (nchar(token) == 0) {
        token <- tryCatch(
            system2("gh", c("auth", "token"), stdout = TRUE, stderr = FALSE),
            error = function(e) ""
        )
    }
    has_token <- nchar(token) > 0
    if (verbose) {
        if (has_token) {
            cli::cli_alert_success("GitHub token available")
        } else {
            cli::cli_alert_warning(
                "No GITHUB_TOKEN found (may hit rate limits during install)"
            )
            cli::cli_text(
                "
    Fix: {.url https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token}"
            )
        }
    }
    has_token
}
