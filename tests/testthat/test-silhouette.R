# Tests the approxSilhouette() function.
# library(testthat); library(bluster); source('test-silhouette.R')

set.seed(80000)
test_that('approxSilhouette yields sensible output for pure clusters', {
    clusters <- rep(seq_len(10), each=100)
    y <- matrix(clusters, ncol=50, nrow=length(clusters))

    out <- approxSilhouette(y, clusters)
    expect_identical(nrow(out), nrow(y))
    expect_true(all(out$width == 1))
    expect_true(all(clusters != out$other))

    # Throwing in some jitter.
    y <- jitter(y)
    out <- approxSilhouette(y, clusters)
    expect_true(all(out$width >= 0.5))
    expect_true(all(clusters != out$other))
})

test_that('approxSilhouette yields correct output for perfectly randomized clusters', {
    clusters <- rep(1:5, each=10)
    y0 <- matrix(rnorm(100), ncol=10)
    y <- rbind(y0, y0, y0, y0, y0)    

    out <- approxSilhouette(y, clusters)
    expect_identical(nrow(out), nrow(y))
    expect_true(all(out$width == 0))
    expect_true(all(clusters != out$other))
})

test_that("approxSilhouette preserves type", {
    y <- matrix(rnorm(1000), ncol=10)
    kclust <- kmeans(y, 5)$cluster
    ref <- approxSilhouette(y, kclust)

    f <- LETTERS[kclust]
    out2 <- approxSilhouette(y, f)
    expect_identical(out2$cluster, f)
    expect_identical(out2$other, LETTERS[ref$other])
    expect_identical(out2$width, ref$width)

    f <- factor(f)
    out2 <- approxSilhouette(y, f)
    expect_identical(out2$cluster, f)
    expect_identical(out2$other, factor(LETTERS[1:5])[ref$other])
    expect_identical(out2$width, ref$width)
})

set.seed(80001)
test_that('approxSilhouette computes the right approximation', {
    y <- matrix(rnorm(1000), ncol=1)
    cout <- clusterRows(y, BLUSPARAM=KmeansParam(4))

    tY <- t(y)
    collated <- numeric(nrow(y))
    for (i in seq_len(nrow(y))) {
        d <- colSums((tY - tY[,i])^2)
        by.clust <- split(d, cout)
        ave.d <- sqrt(vapply(by.clust, mean, 0))

        m <- match(as.character(cout[i]), names(ave.d))
        rest <- ave.d[-m]
        other <- min(rest)
        collated[i] <- (other - ave.d[m])/max(other, ave.d[m])
    } 

    X <- approxSilhouette(y, cout)
    expect_equal(X$width, collated)
    expect_true(all(X$other!=cout))
})
