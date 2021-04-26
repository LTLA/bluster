# Tests the AgnesParam class.
# library(bluster); library(testthat); source('test-agnes-param.R')

test_that("AgnesParam constructor and utilities work correctly", {
    X <- AgnesParam()
    expect_output(show(X), "AgnesParam")

    expect_identical(X[["method"]], NULL)
    X[["method"]] <- "average"
    expect_identical(X[["method"]], "average")

    X <- AgnesParam(cut.params=list(method="average"))
    expect_identical(X[["cut.params"]], list(method="average"))
    X[["cut.params"]] <- list(whee=2)
    expect_identical(X[["cut.params"]], list(whee=2))

    # other show methods
    expect_output(show(AgnesParam(method="average")), "average")
    expect_output(show(AgnesParam(cut.dynamic=TRUE)), "cutreeDynamic")
    expect_output(show(AgnesParam(cut.fun=identity)), "custom")
})

test_that("AgnesParam validity works correctly", {
    expect_error(AgnesParam(cut.dynamic=NA), "non-missing")
    expect_error(AgnesParam(method=1), "character")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, AgnesParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    out2 <- clusterRows(m, AgnesParam(cut.params=list(h=2)))
    expect_identical(length(out2), nrow(m))
    expect_false(identical(out, out2))

    out <- clusterRows(m, AgnesParam(cut.params=list(k=5)))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    ref <- clusterRows(m, AgnesParam(cut.params=list(k=3)))
    out <- clusterRows(m, AgnesParam(cut.fun=function(x) cutree(x, k=3)))
    expect_identical(out, ref)

    full <- clusterRows(m, AgnesParam(), full=TRUE)
    expect_s3_class(full$objects$dist, "dist")
    expect_s3_class(full$objects$hclust, "hclust")

    # Default cut works as expected.
    maxh <- max(full$objects$hclust$height)
    check <- clusterRows(m, AgnesParam(cut.params=list(h=maxh/2)))
    expect_identical(full$clusters, check)
})

test_that("clusterRows works with the dynamic tree cut", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, AgnesParam(cut.dynamic=TRUE))
    expect_true(is.factor(out))
    expect_identical(names(out), NULL)
    expect_identical(length(out), nrow(m))
})

