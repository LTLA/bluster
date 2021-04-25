# Tests the ClaraParam class.
# library(bluster); library(testthat); source('test-clara-param.R')

test_that("ClaraParam constructor and utilities work correctly", {
    X <- ClaraParam(centers=10)
    expect_output(show(X), "ClaraParam")
    expect_output(show(X), "centers: 10")

    expect_identical(X[["centers"]], 10L)
    X[["centers"]] <- 2L
    expect_identical(X[["centers"]], 2L)

    X <- ClaraParam(centers=log)
    expect_true(is.function(X[["centers"]]))
    expect_identical(centers(X, 20), as.integer(round(log(20))))

    X <- ClaraParam(centers=10, sampsize=50)
    expect_identical(X[["sampsize"]], 50)
    X[["sampsize"]] <- 10L
    expect_identical(X[["sampsize"]], 10L)
})

test_that("ClaraParam validity works correctly", {
    expect_error(ClaraParam(-1), "positive")
    expect_error(ClaraParam(10, sampsize=LETTERS), "numeric scalar")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, ClaraParam(5))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    set.seed(9999)
    out <- clusterRows(m, ClaraParam(sqrt))
    expect_identical(length(out), nrow(m))
    expect_equal(nlevels(out), round(sqrt(nrow(m))))

    set.seed(9999)
    full <- clusterRows(m, ClaraParam(sqrt), full=TRUE)
    expect_identical(out, full$cluster)
    expect_s3_class(full$objects$clara, "clara")

    # Responds to the options.
    set.seed(100000)
    suppressWarnings(ref <- cluster::clara(m, k=20, samples=10, sampsize=50, rngR=TRUE, pamLike=TRUE)$clustering)
    set.seed(100000)
    suppressWarnings(out <- clusterRows(m, ClaraParam(centers=20, samples=10, sampsize=50)))
    expect_identical(factor(ref), out)
})
