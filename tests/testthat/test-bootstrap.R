# This tests the bootstrapStability functionality.
# library(testthat); library(bluster); source("test-bootstrap.R")

set.seed(50000)
nobs <- 700
ndims <- 100

set.seed(500001)
test_that("bootstrapStability works correctly with clear separation", {
    dummy <- matrix(rnbinom(nobs*ndims, mu=10, size=20), ncol=ndims)
    known.clusters <- sample(3, nobs, replace=TRUE)
    dummy[known.clusters==1L, 1:30] <- 0
    dummy[known.clusters==2L,31:60] <- 0
    dummy[known.clusters==3L,61:90] <- 0

    output <- bootstrapStability(dummy, FUN=function(x) { kmeans(log10(x+1), 3)$cluster })
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))

    # Works with the mean.
    output <- bootstrapStability(dummy, FUN=function(x) { kmeans(log10(x+1), 3)$cluster }, average="mean")
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))

    # Continues to work if vector is a character or factor.
    output <- bootstrapStability(dummy, FUN=function(x) { c("X", "Y", "Z")[kmeans(log10(x+1), 3)$cluster] })
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))

    output <- bootstrapStability(dummy, FUN=function(x) { factor(kmeans(log10(x+1), 3)$cluster) })
    expect_true(all(output[upper.tri(output, diag=FALSE)] > 0.5))
    expect_true(all(diag(output) > 0.5))
})

set.seed(500002)
test_that("bootstrapStability works correctly with poor separation", {
    dummy <- matrix(rnbinom(nobs*ndims, mu=10, size=20), ncol=ndims)
    output <- bootstrapStability(dummy, FUN=function(x) { kmeans(log10(x+1), 3)$cluster })

    expect_true(all(output[upper.tri(output, diag=TRUE)] < 0.1))
    expect_true(all(output[upper.tri(output, diag=TRUE)] > -0.1))
    expect_true(all(diag(output) < 0.1))
    expect_true(all(diag(output) > -0.1))

    # Works with the mean.
    output <- bootstrapStability(dummy, FUN=function(x) { kmeans(log10(x+1), 3)$cluster }, average="mean")
    expect_true(all(output[upper.tri(output, diag=TRUE)] < 0.1))
    expect_true(all(output[upper.tri(output, diag=TRUE)] > -0.1))
    expect_true(all(diag(output) < 0.1))
    expect_true(all(diag(output) > -0.1))
})

set.seed(500004)
test_that("bootstrapStability works when some clusters are not in the bootstrap.", {
    dummy <- matrix(rnorm(10), nrow=10)

    # Guaranteed to get missing clusters from resampling.
    output <- bootstrapStability(dummy, FUN=function(x) { seq_len(ncol(x)) })

    expect_identical(rownames(output), as.character(seq_len(ncol(dummy))))
    expect_true(all(output[upper.tri(output, diag=FALSE)]==0))
})

set.seed(500003)
test_that("bootstrapStability works with alternative comparison functions", {
    dummy <- matrix(rnorm(nobs*20), nrow=nobs, ncol=20)
    output <- bootstrapStability(dummy, FUN=function(x) { kmeans(x, 3)$cluster }, compare=function(...) 1)
    expect_identical(output, 1)
})

set.seed(500004)
test_that("other miscellaneous tests for bootstrapStability", {
    dummy <- matrix(rnorm(nobs*20), nrow=nobs, ncol=20)

    # Responds correctly to the seed.
    set.seed(20)
    ref <- bootstrapStability(dummy, FUN=function(x) { kmeans(x, 3)$cluster })
    set.seed(20)
    output <- bootstrapStability(dummy, FUN=function(x) { kmeans(x, 3)$cluster })

    expect_identical(ref, output)

    # Errors out.
    expect_error(bootstrapStability(dummy, FUN=function(x) { seq_len(ncol(x)) }, iterations=0), "positive")
})
