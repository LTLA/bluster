# This tests the mergeCommunities function.
# library(testthat); library(bluster); source("test-merge-comm.R")

set.seed(100010)
test_that("mergeCommunities works as expected", {
    m <- matrix(runif(1000), ncol=10)
    output <- clusterRows(m, NNGraphParam(k=5), full=TRUE)
    merged <- mergeCommunities(output$objects$graph, output$clusters, number=3)

    expect_identical(nlevels(merged), 3L)
    expect_identical(length(unique(merged)), 3L)
    expect_true(length(unique(merged)) < length(unique(output$clusters)))

    # Check for nesting.
    freq <- table(paste(merged, output$clusters))
    expect_identical(length(freq), length(unique(output$clusters)))

    merged2 <- mergeCommunities(output$objects$graph, output$clusters, number=2)
    expect_identical(length(unique(merged2)), 2L)
})

set.seed(100011)
test_that("mergeCommunities does, in fact, minimize the modularity", {
    library(igraph)
    groups <- gl(4, 25)
    m <- matrix(runif(1000), ncol=10) + as.numeric(as.integer(groups) %% 2==0 )

    g <- makeSNNGraph(m, k=10)
    ref <- modularity(g, groups, weights=E(g)$weight)

    merged <- mergeCommunities(g, groups, number=2)
    after <- modularity(g, merged, weights=E(g)$weight)
    expect_true(ref < after)

    silly <- modularity(g, gl(2, 50), weights=E(g)$weight)
    expect_true(silly < after)
})

set.seed(100012)
test_that("mergeCommunities behaves with silly inputs", {
    m <- matrix(runif(1000), ncol=10)
    output <- clusterRows(m, NNGraphParam(k=5), full=TRUE)
    expect_error(mergeCommunities(output$objects$graph, output$clusters), "either")
    
    expect_identical(mergeCommunities(output$objects$graph, output$clusters, step=-1), output$clusters)
    expect_identical(mergeCommunities(output$objects$graph, output$clusters, number=100), output$clusters)

    together <- mergeCommunities(output$objects$graph, output$clusters, step=100)
    expect_identical(nlevels(together), 1L)

    expect_identical(mergeCommunities(output$objects$graph, as.integer(output$clusters), step=0), as.integer(output$clusters))
})
