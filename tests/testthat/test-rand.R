# Tests the pairwiseRand() function.
# library(testthat); library(bluster); source('test-rand.R')

m <- matrix(rnorm(1000), ncol=10)

set.seed(8000)
clust1 <- kmeans(m,3)$cluster
clust2 <- kmeans(m,5)$cluster

test_that("pairwiseRand gives the expected output", {
    ratio <- pairwiseRand(clust1, clust2, adjusted=FALSE)

    expect_identical(dim(ratio), c(3L, 3L))
    expect_identical(rownames(ratio), as.character(1:3))
    expect_identical(colnames(ratio), as.character(1:3))

    expect_true(all(is.na(ratio[lower.tri(ratio)])))

    vals <- ratio[upper.tri(ratio, diag=TRUE)]
    expect_true(all(!is.na(vals)))
    expect_true(all(vals >= 0))
    expect_true(all(vals <= 1))

    # Correct values emitted.
    ratio <- pairwiseRand(clust1, clust1, adjusted=FALSE)
    expect_true(all(ratio[upper.tri(ratio, diag=TRUE)]==1))
})

test_that("pairwiseRand mimics the original rand index", {
    full <- pairwiseRand(clust1, clust2, mode="pairs", adjusted=FALSE)

    # Applying a reference calculation.
    status1 <- outer(clust1, clust1, "==")
    status2 <- outer(clust2, clust2, "==")

    keep <- upper.tri(status1)
    status1 <- status1[keep]
    status2 <- status2[keep]

    a <- sum(status1 & status2)
    b <- sum(!status2 & !status1)

    expect_equal(a, sum(diag(full$correct)))
    expect_equal(b, sum(full$correct[upper.tri(full$correct)]))

    rand <- sum(full$correct, na.rm=TRUE)/sum(full$total, na.rm=TRUE)
    expect_identical(rand, (a+b)/choose(length(clust1), 2))
    expect_identical(rand, pairwiseRand(clust1, clust2, mode="index", adjusted=FALSE))
})

test_that("pairwiseRand computes the adjusted Rand index", {
    tab <- table(clust1, clust2)
    same1 <- sum(choose(table(clust1), 2))
    same2 <- sum(choose(table(clust2), 2))
    total <- choose(length(clust1), 2)

    expected <-  same1 * same2 / total
    ref <- (sum(choose(tab, 2)) - expected) / (0.5 * (same1 + same2) - expected)

    ari <- pairwiseRand(clust1, clust2, mode="index", adjusted=TRUE)
    expect_equal(ari, ref)

    # Expected result for equal clustering:
    ari <- pairwiseRand(clust1, clust1, mode="index", adjusted=TRUE)
    expect_equal(ari, 1)
})

test_that("pairwiseRand handles silly inputs", {
    ref <- pairwiseRand(clust1, clust2)

    # Respects empty levels in 'ref'.
    clustX <- factor(clust1, levels=1:5)
    out <- pairwiseRand(clustX, clust2)

    expect_identical(dim(out), c(5L, 5L))
    expect_identical(rownames(out), as.character(1:5))
    expect_identical(colnames(out), as.character(1:5))
    expect_identical(out[1:3,1:3], ref)

    # Ignores empty levels in 'alt'.
    clustY <- factor(clust2, levels=1:5)
    out <- pairwiseRand(clust1, clustY)
    expect_identical(ref, out)

    # Behaves sensibly with zero-length inputs.
    out <- pairwiseRand(clust1[0], clust2[0])
    expect_identical(dim(out), c(0L, 0L))
})

