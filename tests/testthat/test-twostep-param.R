# Tests the TwoStepParam class.
# library(bluster); library(testthat); source('test-twostep-param.R')

test_that("TwoStepParam constructor and utilities work correctly", {
    X <- TwoStepParam()
    expect_output(show(X), "TwoStepParam")

    expect_s4_class(X[["first"]], "KmeansParam")
    X[["first"]] <- HclustParam()
    expect_s4_class(X[["first"]], "HclustParam")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(10000), ncol=10)
    out <- clusterRows(m, TwoStepParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    out2 <- clusterRows(m, TwoStepParam(first=KmeansParam(100)))
    expect_identical(length(out2), nrow(m))
    expect_false(identical(out, out2))

    out <- clusterRows(m, TwoStepParam(second=HclustParam()))
    expect_identical(length(out), nrow(m))
    expect_false(identical(out, out2))

    full <- clusterRows(m, TwoStepParam(), full=TRUE)
    expect_identical(length(full$cluster), nrow(m))
    expect_type(full$objects$first, "list")
    expect_type(full$objects$centroids, "double")
    expect_type(full$objects$second, "list")
})
