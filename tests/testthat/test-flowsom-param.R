# Tests the FlowSOMParam class.
# library(bluster); library(testthat); source('test-flowsom-param.R')

test_that("FlowSOMParam constructor and utilities work correctly", {
    X <- FlowSOMParam(centers=10)
    expect_output(show(X), "FlowSOMParam")
    expect_output(show(X), "centers: 10")
    expect_output(show(X), "alpha: 0.05 0.01")
    expect_output(show(X), "radius: default")
    expect_output(show(X), "rlen: 10")
    expect_output(show(X), "initf: default")
    expect_output(show(X), "distf: 2")

    expect_identical(X[["centers"]], 10L)
    X[["centers"]] <- 2L
    expect_identical(X[["centers"]], 2L)

    X <- FlowSOMParam(centers=log)
    expect_true(is.function(X[["centers"]]))

    X[["radius"]] <- c(0, 0.1)
    expect_output(show(X), "radius: 0 0.1")

    X[["initf"]] <- identity
    expect_output(show(X), "initf: custom")

    X <- FlowSOMParam(centers=10, rlen=50)
    expect_identical(X[["rlen"]], 50L)
    X[["rlen"]] <- 100L 
    expect_identical(X[["rlen"]], 100L)
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(10000), ncol=10)
    colnames(m) <- seq_len(ncol(m))

    set.seed(100)
    out <- clusterRows(m, FlowSOMParam(25))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 25L)

    # Consistent with the defaults in the original function. 
    set.seed(100)
    ref <- FlowSOM::SOM(m, xdim=5, ydim=5)
    expect_identical(out, factor(ref$mapping[,1]))

    # Trying with fewer cells.
    m <- matrix(runif(1000), ncol=10)
    colnames(m) <- seq_len(ncol(m))

    set.seed(100)
    out <- clusterRows(m, FlowSOMParam(12, dim.ratio=3/4))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 12L)

    set.seed(100)
    ref <- FlowSOM::SOM(m, xdim=3, ydim=4)
    expect_identical(out, factor(ref$mapping[,1]))
})

test_that("clusterRows responds to the functions and full=TRUE", {
    m <- matrix(runif(12960), ncol=10)

    set.seed(9999)
    out <- clusterRows(m, FlowSOMParam(sqrt))
    expect_identical(length(out), nrow(m))
    expect_equal(nlevels(out), round(sqrt(nrow(m))))

    set.seed(9999)
    full <- clusterRows(m, FlowSOMParam(sqrt), full=TRUE)
    expect_identical(out, full$cluster)
    expect_true(is.list(full$objects))
})

test_that("clusterRows responds to the options", {
    m <- matrix(runif(10000), ncol=10)
    colnames(m) <- seq_len(ncol(m))

    set.seed(100000)
    suppressWarnings(ref <- FlowSOM::SOM(m, xdim=4, ydim=4, rlen=20, distf=1, init=FALSE, alpha=c(0.12, 0.06)))

    set.seed(100000)
    suppressWarnings(out <- clusterRows(m, FlowSOMParam(16, rlen=20, distf=1, init=FALSE, alpha=c(0.12, 0.06)))) 
    expect_identical(factor(ref$mapping[,1]), out)
})
