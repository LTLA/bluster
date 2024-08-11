# Checks the construction of the SNN graph.
# library(bluster); require(testthat); source("test-make-snn.R")

library(igraph)

check <- function(vals, k=10, type="rank") {
    g <- makeSNNGraph(vals, k=k, type=type) 
    nobs <- nrow(vals)
    expect_identical(seq_len(nobs), as.vector(V(g)))
    nn.out <- BiocNeighbors::findKNN(vals, k=k)

    for (i in seq_len(nobs)) { 
        inn <- c(i, nn.out$index[i,])
        collected <- numeric(nobs)

        for (j in seq_len(nobs)) {
            jnn <- c(j, nn.out$index[j,])
            shared <- intersect(inn, jnn)
            if (length(shared)==0) {
                next
            }
            if (type=="rank") {
                s <- k + 1 - 0.5*(match(shared, inn) + match(shared, jnn))
                collected[j] <- max(c(s, 1e-6))
            } else if (type=="number") {
                collected[j] <- max(length(shared), 1e-6)
            } else {
                collected[j] <- max(length(shared), 1e-6) / length(union(inn, jnn))
            }
        }
        collected[i] <- 0
        expect_equal(collected, g[i])
    }
    return(NULL)
}

set.seed(20000)
nobs <- 200
ndim <- 50

test_that("makeSNNGraph gives same results as a reference", {
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=10)
    
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=20)
    
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=5)

    # Checking 'number' mode.  
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=10, type="number")
    
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=20, type="number")
    
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=5, type="number")

    # Checking 'jaccard' mode.  
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=10, type="jaccard")
    
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=20, type="jaccard")
    
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    check(dummy, k=5, type="jaccard")
})

test_that("makeSNNGraph fails on silly inputs", {
    dummy <- matrix(rnorm(ndim*20), nrow=20)
    expect_warning(out <- makeSNNGraph(dummy, k=50), "capped")
    expect_warning(out2 <- makeSNNGraph(dummy, k=nrow(dummy)-1L), NA)
    expect_identical(out[], out2[])

    expect_error(makeSNNGraph(dummy[,0]), NA) # shouldn't fail, but shouldn't generate anything particularly useful.

    expect_warning(makeSNNGraph(dummy[0,]), "capped")
})

# Checking that makeKNNGraph also works.

KMAKE <- function(dummy, k, directed=FALSE) { 
    nobs <- nrow(dummy)
    collated <- matrix(0, nobs, nobs)
    td <- t(dummy)
    for (cell in seq_len(nobs)) {
        d2 <- colSums((td[,cell] - td)^2)
        chosen <- setdiff(order(d2), cell)[seq_len(k)]
        collated[cell,chosen] <- 1
        if (!directed) { 
            collated[chosen,cell] <- 1
        }
    }
    as(collated, "dgCMatrix")
}

test_that("makeKNNGraph works correctly", {
    dummy <- matrix(rnorm(ndim*nobs), nrow=nobs)
    g <- makeKNNGraph(dummy, k=10)
    expect_false(is.directed(g))
    expect_equal(g[], KMAKE(dummy, k=10))

    g <- makeKNNGraph(dummy, k=10, directed=TRUE)
    expect_true(is.directed(g))
    expect_equal(g[], KMAKE(dummy, k=10, directed=TRUE))

    g <- makeKNNGraph(dummy, k=20)
    expect_equal(g[], KMAKE(dummy, k=20))
})
