# Tests the MbkmeansParam class.
# library(bluster); library(testthat); source('test-mbkmeans-param.R')

test_that("MbkmeansParam constructor and utilities work correctly", {
    X <- MbkmeansParam(centers=10)
    expect_output(show(X), "MbkmeansParam")
    expect_output(show(X), "centers: 10")
    expect_output(show(X), "batch_size: default")
    expect_output(show(X), "BPPARAM: SerialParam")

    expect_identical(X[["centers"]], 10L)
    X[["centers"]] <- 2L
    expect_identical(X[["centers"]], 2L)

    X <- MbkmeansParam(centers=log)
    expect_true(is.function(X[["centers"]]))

    X <- MbkmeansParam(centers=10, batch_size=500)
    expect_identical(X[["batch_size"]], 500L)
    X[["batch_size"]] <- 100L 
    expect_identical(X[["batch_size"]], 100L)
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(10000), ncol=10)

    set.seed(100)
    out <- clusterRows(m, MbkmeansParam(5))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    # Consistent with the defaults in the original function. 
    set.seed(100)
    ref <- mbkmeans::mbkmeans(t(m), 5)
    expect_identical(out, factor(ref$Clusters))

    # Trying with fewer cells.
    m <- matrix(runif(1000), ncol=10)

    set.seed(100)
    out <- clusterRows(m, MbkmeansParam(5))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    set.seed(100)
    ref <- mbkmeans::mbkmeans(t(m), 5)
    expect_identical(out, factor(ref$Clusters))
})

test_that("clusterRows responds to the functions and full=TRUE", {
    m <- matrix(runif(10000), ncol=10)

    set.seed(9999)
    out <- clusterRows(m, MbkmeansParam(sqrt))
    expect_identical(length(out), nrow(m))
    expect_equal(nlevels(out), round(sqrt(nrow(m))))

    set.seed(9999)
    full <- clusterRows(m, MbkmeansParam(sqrt), full=TRUE)
    expect_identical(out, full$cluster)
    expect_true(is.list(full$objects))
})

test_that("clusterRows responds to the options", {
    m <- matrix(runif(10000), ncol=10)
    set.seed(100000)
    suppressWarnings(ref <- mbkmeans::mbkmeans(t(m), 10, batch_size=120, max_iters=5, num_init=10, init_fraction=0.1)$Cluster)

    set.seed(100000)
    suppressWarnings(out <- clusterRows(m, MbkmeansParam(10, batch_size=120, max_iters=5, num_init=10, init_fraction=0.1)))
    expect_identical(factor(ref), out)
})
