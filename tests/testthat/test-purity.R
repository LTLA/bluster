# Tests the neighborPurity() function.
library(testthat); library(bluster); #source('test-purity.R')

set.seed(70000)
test_that('neighborPurity yields correct output for pure clusters', {
    clusters <- rep(seq_len(10), each=100)
    y <- matrix(clusters, ncol=5, nrow=length(clusters))
    y <- jitter(y)

    out <- neighborPurity(y, clusters)
    expect_identical(nrow(out), nrow(y))

    expect_true(all(out$purity==1))
    expect_true(all(out$maximum==clusters))
})

set.seed(70001)
test_that('neighborPurity yields correct output for compromised clusters', {
    y <- matrix(rep(0:1, each=100), ncol=5, nrow=200)
    y <- jitter(y)

    clusters <- rep(1:2, each=100)
    clusters[1] <- 2

    out <- neighborPurity(y, clusters)

    expect_true(out$purity[1] <= 0.05)
    expect_true(all(out$purity[-1] > 0.9))
    expect_identical(out$maximum[1], 1)
    expect_identical(out$maximum[-1], clusters[-1])
})

set.seed(700011)
test_that("neighborPurity handles the weighting correctly", {
    # Creating a bulk of points.
    y0 <- matrix(rnorm(10000), ncol=50)
    y1 <- y2 <- matrix(1000, ncol=50, nrow=10)

    y <- rbind(y1, y2, y0)
    clusters <- rep(1:3, c(nrow(y1), nrow(y2), nrow(y0)))
    out1 <- neighborPurity(y, clusters)$purity
    expect_true(all(abs(out1[1:20]-0.5) < 1e-8))
    expect_true(all(out1[-(1:20)]==1))

    # Unaffected by changes in the number of cells. Technically this should be
    # identical(), but the order of cells returned by findNeighbors() changes,
    # and this results in numeric precision changes due to order of addition of
    # double-precision values, resulting in very slightly different output.
    sub <- c(1:5, nrow(y1) + seq_len(nrow(y2) + nrow(y0)))
    out2 <- neighborPurity(y[sub,], clusters[sub])$purity
    expect_equal(out1[sub], out2)

    sub <- c(1, nrow(y1) + seq_len(nrow(y2) + nrow(y0)))
    out2 <- neighborPurity(y[sub,], clusters[sub])$purity
    expect_equal(out1[sub], out2)
})

set.seed(700012)
test_that("neighborPurity handles other weighting options", {
    # Creating a bulk of points.
    y0 <- matrix(rnorm(10000), ncol=50)
    y1 <- y2 <- matrix(1000, ncol=50, nrow=10)

    y <- rbind(y1, y2, y0)
    clusters <- rep(1:3, c(nrow(y1), nrow(y2), nrow(y0)))

    # Turning off weighting has no effect for balanced clusters.
    out1 <- neighborPurity(y, clusters, weighted=FALSE)$purity
    expect_true(all(abs(out1[1:20]-0.5) < 1e-8))
    expect_true(all(out1[-(1:20)]==1))

    # Turning off weighting has some effect for non-balanced clusters.
    sub <- c(1:5, nrow(y1) + seq_len(nrow(y2) + nrow(y0)))
    out2 <- neighborPurity(y[sub,], clusters[sub], weighted=FALSE)$purity
    expect_true(all(abs(out2[1:5]-1/3) < 1e-8))
    expect_true(all(abs(out2[6:15]-2/3) < 1e-8))
    expect_true(all(out2[-(1:20)]==1))

    # We can replace it with our own weighting to restore the balance.
    out3 <- neighborPurity(y[sub,], clusters[sub], weighted=rep(c(2, 1), c(5, nrow(y2)+nrow(y0))))$purity
    ref <- neighborPurity(y[sub,], clusters[sub])$purity
    expect_equal(out3, ref)
})
