# Tests the DianaParam class.
# library(bluster); library(testthat); source('test-diana-param.R')

test_that("DianaParam constructor and utilities work correctly", {
    X <- DianaParam()
    expect_output(show(X), "DianaParam")

    expect_identical(X[["metric"]], NULL)
    X[["metric"]] <- "euclidean"
    expect_identical(X[["metric"]], "euclidean")

    X <- DianaParam(cut.params=list(metric="euclidean"))
    expect_identical(X[["cut.params"]], list(metric="euclidean"))
    X[["cut.params"]] <- list(whee=2)
    expect_identical(X[["cut.params"]], list(whee=2))

    # other show metrics
    expect_output(show(DianaParam(metric="euclidean")), "euclidean")
    expect_output(show(DianaParam(cut.dynamic=TRUE)), "cutreeDynamic")
    expect_output(show(DianaParam(cut.fun=identity)), "custom")
})

test_that("DianaParam validity works correctly", {
    expect_error(DianaParam(cut.dynamic=NA), "non-missing")
    expect_error(DianaParam(metric=1), "character")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, DianaParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    out2 <- clusterRows(m, DianaParam(cut.params=list(h=2)))
    expect_identical(length(out2), nrow(m))
    expect_false(identical(out, out2))

    out <- clusterRows(m, DianaParam(cut.params=list(k=5)))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    ref <- clusterRows(m, DianaParam(cut.params=list(k=3)))
    out <- clusterRows(m, DianaParam(cut.fun=function(x) cutree(x, k=3)))
    expect_identical(out, ref)

    full <- clusterRows(m, DianaParam(), full=TRUE)
    expect_s3_class(full$objects$dist, "dist")
    expect_s3_class(full$objects$hclust, "hclust")

    # Default cut works as expected.
    maxh <- max(full$objects$hclust$height)
    check <- clusterRows(m, DianaParam(cut.params=list(h=maxh/2)))
    expect_identical(full$clusters, check)
})

test_that("clusterRows works with the dynamic tree cut", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, DianaParam(cut.dynamic=TRUE))
    expect_true(is.factor(out))
    expect_identical(names(out), NULL)
    expect_identical(length(out), nrow(m))
})

test_that("clusterRows works with custom distance functions", {
    m <- matrix(runif(1000), ncol=10)

    out <- clusterRows(m, DianaParam(metric="manhattan"), full=TRUE)
    expect_equal(as.matrix(out$objects$dist), as.matrix(dist(m, method="manhattan")))

    vegout <- clusterRows(m, DianaParam(metric = "canberra", dist.fun = vegan::vegdist), full=TRUE)
    expected <- vegan::vegdist(m, "canberra")
    expect_equal(as.matrix(vegout$objects$dist), as.matrix(expected))
})

