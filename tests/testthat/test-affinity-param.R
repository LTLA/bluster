# Tests the AffinityParam class.
# library(bluster); library(testthat); source('test-affinity-param.R')

set.seed(1000)
test_that("AffinityParam constructor and utilities work correctly", {
    X <- AffinityParam()
    expect_output(show(X), "AffinityParam")
    expect_output(show(X), "s: default")
    expect_output(show(X), "p: NA")

    expect_identical(X[["s"]], NULL)
    X[["s"]] <- identity
    expect_identical(X[["s"]], identity)

    X <- AffinityParam(q=0.9)
    expect_identical(X[["q"]], 0.9)
    X[["q"]] <- 0.1
    expect_identical(X[["q"]], 0.1)
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(1000), ncol=10)

    out <- clusterRows(m, AffinityParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    # Consistent with reference values.
    ref <- apcluster::apcluster(apcluster::negDistMat(m, r=2))
    expect_identical(out, factor(apcluster::labels(ref, type="enum")))

    # Trying with different parameters.
    set.seed(10)
    out <- clusterRows(m, AffinityParam(q=0.8))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))

    set.seed(10)
    ref <- apcluster::apcluster(apcluster::negDistMat(m, r=2), q=0.8)
    expect_identical(out, factor(apcluster::labels(ref, type="enum")))
})

test_that("clusterRows responds to full=TRUE", {
    m <- matrix(runif(2960), ncol=10)
    out <- clusterRows(m, AffinityParam())

    full <- clusterRows(m, AffinityParam(), full=TRUE)
    expect_identical(out, full$cluster)
    expect_s4_class(full$objects$apcluster, "APResult")
})

set.seed(1010001)
test_that("clusterRows responds to the options", {
    m <- matrix(runif(1000), ncol=10)

    set.seed(1000)
    out <- clusterRows(m, AffinityParam(q=0.8, lam=0.7, s=apcluster::expSimMat()))

    set.seed(1000)
    ref <- apcluster::apcluster(apcluster::expSimMat(m), q=0.8, lam=0.7)
    expect_identical(out, factor(apcluster::labels(ref, type="enum")))

    # Responds to the negative q.
    set.seed(2000)
    out <- clusterRows(m, AffinityParam(q=-1, s=apcluster::expSimMat()))

    set.seed(2000)
    mat <- apcluster::expSimMat(m)
    ref <- apcluster::apcluster(mat, p=0)
    expect_identical(out, factor(apcluster::labels(ref, type="enum")))

    # Also works for negative distances.
    set.seed(300)
    out <- clusterRows(m, AffinityParam(q=-2, s=apcluster::negDistMat()))

    set.seed(300)
    mat <- apcluster::negDistMat(m)
    ref <- apcluster::apcluster(mat, p=min(mat) * 3)
    expect_identical(out, factor(apcluster::labels(ref, type="enum")))
})
