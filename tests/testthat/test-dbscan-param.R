# Tests the DbscanParam class.
# library(bluster); library(testthat); source('test-dbscan-param.R')

set.seed(1000)
test_that("DbscanParam constructor and utilities work correctly", {
    X <- DbscanParam()
    expect_output(show(X), "DbscanParam")
    expect_output(show(X), "eps: default")
    expect_output(show(X), "min.pts: 5")

    expect_identical(X[["eps"]], NULL)
    X[["eps"]] <- 5 
    expect_identical(X[["eps"]], 5)

    X <- DbscanParam(min.pts=10)
    expect_identical(X[["min.pts"]], 10L)
    X[["min.pts"]] <- 2L
    expect_identical(X[["min.pts"]], 2L)
})

set.seed(10)
raw <- matrix(rnorm(100), ncol=10)
m <- matrix(raw[sample(nrow(raw), 1000, replace=TRUE),], ncol=10)
m <- m + rnorm(length(m), sd=0.1)

# Creating an alternative implementation with less subsetting all over the place.
REF <- function(m, eps, min.pts) {
    # Identifying the core points.
    neighbors <- BiocNeighbors::findNeighbors(m, threshold=eps, get.distance=FALSE)$index
    is.core <- lengths(neighbors) > min.pts # >, not >=, as self is included in range.
    core.id <- which(is.core)
    core.neighbors <- neighbors[is.core]

    # Constructing the core components.
    from <- rep(core.id, lengths(core.neighbors))
    to <- unlist(core.neighbors)
    keep <- to %in% from
    g <- igraph::make_graph(rbind(as.character(from[keep]), as.character(to[keep])), directed=FALSE)

    comp <- igraph::components(g)$membership
    clusters <- integer(nrow(m))
    clusters[as.integer(names(comp))] <- unname(comp)

    # Matching all points to their closest core point.
    closest.core <- BiocNeighbors::queryKNN(X=m[is.core,,drop=FALSE], query=m, k=1)
    reassign <- closest.core$distance <= eps
    clusters[reassign] <- comp[as.character(core.id)][closest.core$index[reassign]]

    factor(ifelse(clusters == 0, NA_integer_, clusters))
}

test_that("clusterRows works correctly", {
    out <- clusterRows(m, DbscanParam())
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 10L)

    # Trying with different parameters.
    out2 <- clusterRows(m, DbscanParam(core.prop=0.8))
    expect_true(is.factor(out2))
    expect_identical(length(out2), nrow(m))
    expect_true(sum(is.na(out2)) < sum(is.na(out)))

    # Comparing to another implementation.
    out3 <- clusterRows(m, DbscanParam(eps=0.3, min.pts=5))
    expect_identical(nlevels(out3), 10L)
    ref <- REF(m, eps=0.3, min.pts=5)
    expect_identical(out3, ref)

    out4 <- clusterRows(m, DbscanParam(eps=0.5, min.pts=10))
    expect_identical(nlevels(out4), 10L)
    ref <- REF(m, eps=0.5, min.pts=10)
    expect_identical(out4, ref)

    out5 <- clusterRows(m, DbscanParam(eps=0.25, min.pts=10))
    expect_true(nlevels(out5) > 5)
    ref <- REF(m, eps=0.25, min.pts=10)
    expect_identical(out5, ref)

    # Edge case behaves well.
    none <- clusterRows(m, DbscanParam(eps=0.001))
    expect_true(all(is.na(none)))
})

test_that("clusterRows responds to full=TRUE", {
    out <- clusterRows(m, DbscanParam())
    full <- clusterRows(m, DbscanParam(), full=TRUE)
    expect_identical(out, full$cluster)
    expect_type(full$objects$min.pts, "integer")
    expect_type(full$objects$eps, "double")
})
