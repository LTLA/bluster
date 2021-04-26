# This tests the sweeping functionality.
# library(testthat); library(bluster); source("test-sweep.R")

test_that("clusterSweep works correctly in simple cases", {
    x <- matrix(rnorm(1000), ncol=5)
    out <- clusterSweep(x, NNGraphParam(), k=c(5L, 10L, 20L), cluster.fun=c("louvain", "walktrap"))
    for (choice in seq_len(nrow(out$parameters))) {
        curk <- out$parameters$k[[choice]]
        curf <- out$parameters$cluster.fun[[choice]]
        param <- NNGraphParam(k=curk, cluster.fun=curf)
        expect_identical(out$clusters[,choice], clusterRows(x, param))
        expect_identical(colnames(out$clusters)[choice], sprintf("k.%s_cluster.fun.%s", curk, curf))
    }

    # Respects existing parameters.
    x <- matrix(rnorm(1000), ncol=5)
    out <- clusterSweep(x, NNGraphParam(k=12), cluster.fun=c("louvain", "walktrap"))
    for (choice in 1:2) {
        curf <- out$parameters$cluster.fun[[choice]]
        param <- NNGraphParam(k=12, cluster.fun=curf)
        expect_identical(out$clusters[,choice], clusterRows(x, param))
        expect_identical(colnames(out$clusters)[choice], sprintf("cluster.fun.%s", curf))
    }
})

test_that("clusterSweep reports full objects correctly", {
    x <- matrix(rnorm(1000), ncol=5)

    out <- clusterSweep(x, NNGraphParam(), cluster.fun=c("louvain", "walktrap"), full=TRUE)
    ref <- clusterSweep(x, NNGraphParam(), cluster.fun=c("louvain", "walktrap"))
    expect_identical(out[1:2], ref)

    counter <- clusterRows(x, NNGraphParam(cluster.fun="walktrap"), full=TRUE)
    expect_identical(out$objects$cluster.fun.walktrap$communities, counter$objects$communities)
})

test_that("clusterSweep handles difficult names correctly", {
    x <- matrix(rnorm(1000), ncol=5)
    out <- clusterSweep(x, NNGraphParam(cluster.fun="walktrap"), cluster.args=list(list(steps=3), list(steps=4)))
    expect_identical(colnames(out$clusters), rep("cluster.args.list", 2))
})

test_that("clusterSweep loops across list elements correctly", {
    set.seed(100)
    x <- matrix(rnorm(1000), ncol=5)
    out <- clusterSweep(x, NNGraphParam(cluster.fun="walktrap"), cluster.args=list(list(steps=3), list(steps=4)))

    # Checking that the list is actually extracted properly.
    ref <- clusterRows(x, NNGraphParam(cluster.fun="walktrap"))
    expect_identical(out$clusters[,2], ref)
    expect_false(identical(out$clusters[,1], ref)) # actually has an effect.

    alt <- clusterSweep(x, NNGraphParam(cluster.fun="walktrap"), cluster.args=list(steps=3:4))
    expect_equal(out, alt)

    expanded <- clusterSweep(x[1:20,], NNGraphParam(cluster.fun="spinglass", k=2), cluster.args=list(spins=c(20, 25, 30), start.temp=c(0.9, 1, 1.1)))
    expect_identical(nrow(expanded$parameters), 9L)
    expect_identical(expanded$parameters[1,1][[1]], list(spins=20, start.temp=0.9))
    expect_identical(expanded$parameters[9,1][[1]], list(spins=30, start.temp=1.1))
})

test_that("clusterSweep skips looping over AsIs'd list elements", {
    out <- bluster:::.expand_list(list(k=1:50))
    expect_identical(length(out), 50L)

    out <- bluster:::.expand_list(list(k=I(1:50)))
    expect_identical(out, list(list(k=1:50)))

    out <- bluster:::.expand_list(list(X=HclustParam()))
    expect_identical(out, list(list(X=HclustParam())))
})
