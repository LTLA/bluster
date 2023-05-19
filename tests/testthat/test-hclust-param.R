# Tests the HclustParam class.
# library(bluster); library(testthat); source('test-hclust-param.R')

test_that("HclustParam constructor and utilities work correctly", {
    X <- HclustParam()
    expect_output(show(X), "HclustParam")

    expect_identical(X[["method"]], NULL)
    X[["method"]] <- "average"
    expect_identical(X[["method"]], "average")

    X <- HclustParam(cut.params=list(method="average"))
    expect_identical(X[["cut.params"]], list(method="average"))
    X[["cut.params"]] <- list(whee=2)
    expect_identical(X[["cut.params"]], list(whee=2))

    # other show methods
    expect_output(show(HclustParam(method="average")), "average")
    expect_output(show(HclustParam(cut.dynamic=TRUE)), "cutreeDynamic")
    expect_output(show(HclustParam(cut.fun=identity)), "custom")
})

test_that("HclustParam validity works correctly", {
    expect_error(HclustParam(cut.dynamic=NA), "non-missing")
    expect_error(HclustParam(method=1), "character")
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, HclustParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    out2 <- clusterRows(m, HclustParam(cut.params=list(h=2)))
    expect_identical(length(out2), nrow(m))
    expect_false(identical(out, out2))

    out <- clusterRows(m, HclustParam(cut.params=list(k=5)))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 5L)

    ref <- clusterRows(m, HclustParam(cut.params=list(k=3)))
    out <- clusterRows(m, HclustParam(cut.fun=function(x) cutree(x, k=3)))
    expect_identical(out, ref)

    full <- clusterRows(m, HclustParam(), full=TRUE)
    expect_s3_class(full$objects$dist, "dist")
    expect_s3_class(full$objects$hclust, "hclust")

    # Default cut works as expected.
    maxh <- max(full$objects$hclust$height)
    check <- clusterRows(m, HclustParam(cut.params=list(h=maxh/2)))
    expect_identical(full$clusters, check)
})

test_that("clusterRows works with the dynamic tree cut", {
    m <- matrix(runif(1000), ncol=10)
    out <- clusterRows(m, HclustParam(cut.dynamic=TRUE))
    expect_true(is.factor(out))
    expect_identical(names(out), NULL)
    expect_identical(length(out), nrow(m))
})

# The following test is only ran when vegan is available
if (require(vegan)) {
    test_that("dist.fun parameter works correctly", {
        m <- matrix(runif(1000), ncol=10)
        dist_result <- clusterRows(m, HclustParam(metric = "euclidean"))
        vegdist_result <- clusterRows(m, HclustParam(metric = "euclidean", dist.fun = vegan::vegdist))
        expect_identical(dist_result, vegdist_result)
        
        vegdist_result2 <- clusterRows(m, HclustParam(metric = "canberra", dist.fun = vegan::vegdist), full = TRUE)
        res_dist_matrix <- vegdist_result2$objects$dist
        original_dist_matrix <- vegan::vegdist(m, "canberra")
        expect_setequal(res_dist_matrix, original_dist_matrix)
    })
}