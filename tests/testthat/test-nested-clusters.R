# This contains some cursory tests for nestedClusters().
# library(testthat); library(bluster); source("test-nested-clusters.R")

set.seed(1000)
m <- matrix(runif(10000), ncol=10)
clust1 <- kmeans(m,10)$cluster
clust2 <- kmeans(m,20)$cluster

test_that("nestedClusters works as expected", {
    out <- nestedClusters(clust1, clust2)
    expect_identical(rownames(out$alt.mapping), as.character(sort(unique(clust2))))
    expect_identical(out$alt.mapping$max, unname(apply(out$proportions, 1, max)))
    expect_identical(out$alt.mapping$which, as.character(apply(out$proportions, 1, which.max)))

    # Works with some more interesting names.
    out <- nestedClusters(letters[clust1], LETTERS[clust2])
    expect_identical(rownames(out$alt.mapping), LETTERS[1:20])
    expect_identical(out$alt.mapping$which, letters[apply(out$proportions, 1, which.max)])

    # Works with factors.
    out <- nestedClusters(factor(clust1, 10:1), factor(clust2, 20:1))
    expect_identical(rownames(out$alt.mapping), as.character(20:1))
    expect_identical(out$alt.mapping$which, as.character(10:1)[apply(out$proportions, 1, which.max)])
})

test_that("nestedClusters ref.score behaves correctly", {
    # The ref.score is 1 in cases of perfect nesting.
    out <- nestedClusters(clust1, clust1)$ref.score
    expect_true(all(out==1))
                                                                                
    nest.clust <- paste0(clust1, sample(letters, length(clust1), replace=TRUE))
    out <- nestedClusters(clust1, nest.clust)$ref.score
    expect_true(all(out==1))
                                                                            
    # In contrast, it is much lower when nesting is bad.
    out <- nestedClusters(clust1, sample(clust1))$ref.score
    expect_true(all(out < 0.2))
})
