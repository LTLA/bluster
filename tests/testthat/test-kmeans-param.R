# Tests the KmeansParam class.
# library(bluster); library(testthat); source('test-kmeans-param.R')

test_that("KmeansParam constructor and utilities work correctly", {
    X <- KmeansParam(centers=10)
    expect_output(show(X), "KmeansParam")
    expect_output(show(X), "centers: 10")

    expect_identical(X[["centers"]], 10L)
    X[["centers"]] <- 2L
    expect_identical(X[["centers"]], 2L)

    X <- KmeansParam(centers=log)
    expect_true(is.function(X[["centers"]]))

    X <- KmeansParam(centers=10, algorithm="Lloyd")
    expect_identical(X[["algorithm"]], "Lloyd")
    X[["algorithm"]] <- "Hartigan"
    expect_identical(X[["algorithm"]], "Hartigan")
})

test_that("KmeansParam validity works correctly", {
    expect_error(KmeansParam(-1), "positive")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, KmeansParam(5))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    set.seed(9999)
    out <- clusterRows(m, KmeansParam(sqrt))
    expect_identical(length(out), nrow(m))
    expect_equal(nlevels(out), round(sqrt(nrow(m))))

    set.seed(9999)
    full <- clusterRows(m, KmeansParam(sqrt), full=TRUE)
    expect_identical(out, full$cluster)
    expect_s3_class(full$objects, "kmeans")

    # Responds to the options.
    set.seed(100000)
    suppressWarnings(ref <- kmeans(m, centers=20, nstart=5, iter.max=3, algorithm="Lloyd")$cluster)
    set.seed(100000)
    suppressWarnings(out <- clusterRows(m, KmeansParam(centers=20, nstart=5, iter.max=3, algorithm="Lloyd")))
    expect_identical(factor(ref), out)
})
