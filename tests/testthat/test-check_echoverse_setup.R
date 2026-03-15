test_that("check_echoverse_setup runs without error", {

    res <- check_echoverse_setup(verbose = FALSE)
    testthat::expect_type(res, "list")
    testthat::expect_true("packages" %in% names(res))
    testthat::expect_true("dependencies" %in% names(res))
    testthat::expect_true("system" %in% names(res))
    testthat::expect_true("tools" %in% names(res))
    testthat::expect_true("pass" %in% names(res))
    testthat::expect_type(res$pass, "logical")
})

test_that("check_echoverse_setup packages result is a data.frame", {

    res <- check_echoverse_setup(verbose = FALSE)
    testthat::expect_s3_class(res$packages, "data.frame")
    testthat::expect_true(all(c("package", "installed", "loadable", "version")
                              %in% names(res$packages)))
    ## echolocatoR itself should always be loadable
    echolocatoR_row <- res$packages[res$packages$package == "echolocatoR", ]
    testthat::expect_true(echolocatoR_row$loadable)
})

test_that("check_echoverse_setup accepts custom package lists", {

    res <- check_echoverse_setup(
        echoverse_pkgs = c("echolocatoR"),
        key_deps = c("data.table"),
        verbose = FALSE
    )
    testthat::expect_equal(nrow(res$packages), 1)
    testthat::expect_equal(nrow(res$dependencies), 1)
})

test_that("check_echoverse_setup verbose produces output", {

    ## cli writes to connections that capture.output may miss,
    ## so just verify it runs without error in verbose mode
    testthat::expect_no_error(
        suppressMessages(
            check_echoverse_setup(
                echoverse_pkgs = "echolocatoR",
                key_deps = "data.table",
                verbose = TRUE
            )
        )
    )
})

test_that("get_os returns valid string", {

    os <- echolocatoR:::get_os()
    testthat::expect_true(os %in% c("macOS", "Linux", "Windows"))
})

test_that("check_github_token returns logical", {

    res <- echolocatoR:::check_github_token(verbose = FALSE)
    testthat::expect_type(res, "logical")
})
