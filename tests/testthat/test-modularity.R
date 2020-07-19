# Tests the pairwiseModularity function.
# library(bluster); library(testthat); source("test-modularity.R")

set.seed(20001)
test_that("pairwiseModularity computes the correct values", {
    exprs <- matrix(rnorm(1000), ncol=10)
    g <- makeSNNGraph(exprs)

    random <- sample(5, nrow(exprs), replace=TRUE)
    out <- pairwiseModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random, weight=igraph::E(g)$weight))
    expect_equal(sum(out), 0)

    # Some basic checks on the expected values.
    out <- pairwiseModularity(g, random, get.weights=TRUE)
    expect_equal(sum(out$observed), sum(out$expected))
    expect_equal(sum(out$observed), sum(igraph::E(g)$weight))
})

set.seed(20001)
test_that("pairwiseModularity computes correct values in a more realistic scenario", {
    exprs <- matrix(rnorm(5000), ncol=50)
    exprs[ 1:20, 1:10] <- rnorm(200, 2)
    exprs[21:50,11:20] <- rnorm(300, 2)
    exprs[51:90,21:30] <- rnorm(400, 2)
    g <- makeSNNGraph(exprs, k=10)

    actual <- igraph::cluster_walktrap(g)
    out <- pairwiseModularity(g, actual$membership) 
    expect_equal(sum(diag(out)), igraph::modularity(g, actual$membership, weight=igraph::E(g)$weight))
    expect_equal(sum(out), 0)

    # Some basic checks on the expected values.
    out <- pairwiseModularity(g, actual$membership, get.weights=TRUE)
    expect_equal(sum(out$observed), sum(out$expected))
    expect_equal(sum(out$observed), sum(igraph::E(g)$weight))
})

set.seed(20002)
test_that("pairwiseModularity handles unweighted graphs correctly", {
    exprs <- matrix(rnorm(1000), ncol=10)
    g <- makeKNNGraph(exprs)

    random <- sample(5, nrow(exprs), replace=TRUE)
    out <- pairwiseModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random))
    expect_equal(sum(out), 0)

    # Some basic checks on the expected values.
    out <- pairwiseModularity(g, random, get.weights=TRUE)
    expect_equal(sum(out$observed), sum(out$expected))
    expect_equal(sum(out$observed), length(igraph::E(g)))
})

set.seed(20003)
test_that("pairwiseModularity handles directed graphs correctly", {
    exprs <- matrix(rnorm(1000), ncol=10)
    g <- makeKNNGraph(exprs, directed=TRUE)

    random <- sample(5, nrow(exprs), replace=TRUE)
    out <- pairwiseModularity(g, random)
    ref <- pairwiseModularity(igraph::as.undirected(g, mode="each"), random) 
    expect_identical(out, ref)
})

set.seed(20003)
test_that("pairwiseModularity handles self-loops correctly", {
    exprs <- matrix(rnorm(1000), ncol=10)
    g <- makeSNNGraph(exprs)

    g <- igraph::add_edges(g, rep(nrow(exprs), each=2), weight=10)
    random <- sample(5, nrow(exprs), replace=TRUE)
    out <- pairwiseModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random, weight=igraph::E(g)$weight))

    # Works for unweighted graphs.
    exprs <- matrix(rnorm(10000), ncol=10)
    g <- makeKNNGraph(exprs)

    g <- igraph::add_edges(g, rep(nrow(exprs), each=2))
    random <- sample(5, nrow(exprs), replace=TRUE)
    out <- pairwiseModularity(g, random) 
    expect_equal(sum(diag(out)), igraph::modularity(g, random))
})
