# This tests the clusterRMSD function.
# library(testthat); library(bluster); source("test-cluster-rmsd.R")

set.seed(1000)
x <- matrix(rnorm(10000), ncol=10)
kout <- kmeans(x, 5)

test_that("clusterRMSD works as expected", {
    expect_equal(kout$withinss, unname(clusterRMSD(x, kout$cluster, sum=TRUE)))
    expect_identical(names(clusterRMSD(x, kout$cluster)), as.character(1:5))

    ref <- kout$withinss / as.integer(table(kout$cluster) - 1)
    expect_equal(ref, unname(clusterRMSD(x, kout$cluster)))

    out <- clusterRMSD(x, letters[kout$cluster])
    expect_equal(ref, unname(out))
    expect_identical(names(out), letters[1:5])

    expect_identical(clusterRMSD(x[0,], numeric(0)), numeric(0))
})
