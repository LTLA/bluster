# Tests the DMMParam class.
# library(bluster); library(testthat); source('test-dmm-param.R')

set.seed(1000)
test_that("DMMParam constructor and utilities work correctly", {
    X <- DMMParam()
    expect_output(show(X), "DMMParam")
    expect_output(show(X), "k: 1:3")
    expect_output(show(X), "type: laplace")
    expect_output(show(X), "transposed: FALSE")
    
    expect_identical(X[["k"]], 1:3)
    X[["k"]] <- as.integer(2:4)
    expect_identical(X[["k"]], 2:4)
    
    X <- DMMParam(k=3)
    expect_equal(X[["k"]], 3)
    X[["k"]] <- as.integer(2)
    expect_equal(X[["k"]], 2)
    
    X <- DMMParam(type="BIC")
    expect_identical(X[["type"]], "BIC")
    X[["type"]] <- "AIC"
    expect_identical(X[["type"]], "AIC")
})


test_that("DMMParam validity works correctly", {
    expect_error(DMMParam(type="unknown"), "must be equal")
    expect_error(DMMParam(k=-1), "strictly positive")
    expect_error(DMMParam(k=0:3), "strictly positive")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(10000), ncol = 10, nrow = 1000)
    m <- ceiling(m * 10)
    
    out <- clusterRows(m, DMMParam(), full = TRUE)
    expect_true(is.factor(out))
    expect_identical(length(out), ncol(m))
    
    # Consistent with reference values.
    k=out$objects$k
    dmm <- .get_dmm(t(m), k=k)
    prob <- DirichletMultinomial::mixture(dmm[[1]])
    colnames(prob) <- 1:k
    clusters <- colnames(prob)[max.col(prob, ties.method = "first")]
    clusters <- factor(clusters)
    names(clusters) <- rownames(x)
    expect_identical(out, clusters)
    
    # Trying with different parameters.
    set.seed(10)
    out <- clusterRows(m, DMMParam(k=4))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    
    set.seed(10)
    dmm <- .get_dmm(t(m), k=4)
    prob <- DirichletMultinomial::mixture(dmm[[1]])
    colnames(prob) <- 1:4
    clusters <- colnames(prob)[max.col(prob, ties.method = "first")]
    clusters <- factor(clusters)
    names(clusters) <- rownames(x)
    expect_identical(out, clusters)
})

test_that("clusterRows responds to full=TRUE", {
    m <- matrix(runif(1000), ncol = 100, nrow = 10)
    m <- ceiling(m * 10)
    out <- clusterRows(m, DMMParam())
    
    full <- clusterRows(m, DMMParam(), full=TRUE)
    expect_identical(out, full$clusters)
    expect_s4_class(full$objects$dmm[[1]], "DMN")
})
