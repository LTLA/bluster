# Tests the HclustParam class.
# library(bluster); library(testthat); source('test-hclust-param.R')

test_that("HclustParam constructor and utilities work correctly", {
    X <- HclustParam()
    expect_output(show(X), "HclustParam")

    expect_identical(X[["method"]], "complete")
    X[["method"]] <- "average"
    expect_identical(X[["method"]], "average")

    X <- HclustParam(BLAH=2)
    expect_identical(X[["BLAH"]], 2)
    X[["BLAH"]] <- "average"
    expect_identical(X[["BLAH"]], "average")

    # other show methods
    expect_output(show(HclustParam(cut.number=2)), "cut.number")
    expect_output(show(HclustParam(cut.fun=identity)), "custom")
})

test_that("HclustParam validity works correctly", {
    expect_error(HclustParam(NA_character_), "non-missing")
    expect_error(HclustParam(cut.height=-1), "positive")
    expect_error(HclustParam(cut.number=-1), "positive")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, HclustParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    out2 <- clusterRows(m, HclustParam(cut.height=2))
    expect_identical(length(out2), nrow(m))
    expect_false(identical(out, out2))

    out <- clusterRows(m, HclustParam(cut.number=5))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    ref <- clusterRows(m, HclustParam(cut.number=3))
    out <- clusterRows(m, HclustParam(cut.fun=function(x) cutree(x, k=3)))
    expect_identical(out, ref)

    full <- clusterRows(m, HclustParam(cut.number=3), full=TRUE)
    expect_identical(ref, full$cluster)
    expect_s3_class(full$objects, "hclust")
})
