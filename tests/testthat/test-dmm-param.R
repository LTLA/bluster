# Tests the DmmParam class.
# library(bluster); library(testthat); source('test-dmm-param.R')

set.seed(1000)
test_that("DmmParam constructor and utilities work correctly", {
    X <- DmmParam()
    expect_output(show(X), "DmmParam")
    expect_output(show(X), "k(3): 1 2 3", fixed=TRUE)
    expect_output(show(X), "type: laplace")

    expect_identical(X[["k"]], 1:3)
    X[["k"]] <- as.integer(2:4)
    expect_identical(X[["k"]], 2:4)

    X <- DmmParam(k=3)
    expect_equal(X[["k"]], 3)
    X[["k"]] <- as.integer(2)
    expect_equal(X[["k"]], 2)

    X <- DmmParam(type="BIC")
    expect_identical(X[["type"]], "BIC")
    X[["type"]] <- "AIC"
    expect_identical(X[["type"]], "AIC")
})


test_that("DmmParam validity works correctly", {
    expect_error(DmmParam(type="unknown"), "should be one of")
    expect_error(DmmParam(k=-1), "strictly positive")
    expect_error(DmmParam(k=0:3), "strictly positive")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol = 100, nrow = 10)
    m <- ceiling(m * 10)
    rownames(m) <- paste0("A",1:10)
    colnames(m) <- 1:100

    out <- clusterRows(m, DmmParam(), full = TRUE)
    expect_true(is.factor(out$clusters))
    expect_identical(length(out$clusters), nrow(m))

    # Trying with different parameters.
    set.seed(10)
    out <- clusterRows(m, DmmParam(k=4))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
})

test_that("clusterRows responds to full=TRUE", {
    m <- matrix(runif(1000), ncol = 100, nrow = 10)
    m <- ceiling(m * 10)
    rownames(m) <- paste0("A",1:10)
    colnames(m) <- 1:100
    out <- clusterRows(m, DmmParam())

    full <- clusterRows(m, DmmParam(), full=TRUE)
    expect_identical(out, full$clusters)
    expect_s4_class(full$objects$dmm[[1]], "DMN")
})
