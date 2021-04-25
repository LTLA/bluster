# Tests the PamParam class.
# library(bluster); library(testthat); source('test-pam-param.R')

test_that("PamParam constructor and utilities work correctly", {
    X <- PamParam(centers=10)
    expect_output(show(X), "PamParam")
    expect_output(show(X), "centers: 10")

    expect_identical(X[["centers"]], 10L)
    X[["centers"]] <- 2L
    expect_identical(X[["centers"]], 2L)

    X <- PamParam(centers=log)
    expect_true(is.function(X[["centers"]]))
    expect_identical(centers(X, 20), as.integer(round(log(20))))

    X <- PamParam(centers=10, variant="faster")
    expect_identical(X[["variant"]], "faster")
    X[["variant"]] <- "original"
    expect_identical(X[["variant"]], "original")
})

test_that("PamParam validity works correctly", {
    expect_error(PamParam(-1), "positive")
    expect_error(PamParam(10, variant=LETTERS), "character scalar")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, PamParam(5))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    set.seed(9999)
    out <- clusterRows(m, PamParam(sqrt))
    expect_identical(length(out), nrow(m))
    expect_equal(nlevels(out), round(sqrt(nrow(m))))

    set.seed(9999)
    full <- clusterRows(m, PamParam(sqrt), full=TRUE)
    expect_identical(out, full$cluster)
    expect_s3_class(full$objects$pam, "pam")

    # Responds to the options.
    set.seed(100000)
    suppressWarnings(ref <- cluster::pam(m, k=20, variant="faster", do.swap=FALSE)$clustering)
    set.seed(100000)
    suppressWarnings(out <- clusterRows(m, PamParam(centers=20, variant="faster", do.swap=FALSE)))
    expect_identical(factor(ref), out)
})
