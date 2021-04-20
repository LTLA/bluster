# Tests the SOMParam class.
# library(bluster); library(testthat); source('test-som-param.R')

set.seed(1000)
test_that("SOMParam constructor and utilities work correctly", {
    X <- SOMParam(centers=10)
    expect_output(show(X), "SOMParam")
    expect_output(show(X), "centers: 10")
    expect_output(show(X), "alpha: 0.05 0.01")
    expect_output(show(X), "radius: default")
    expect_output(show(X), "rlen: 10")
    expect_output(show(X), "dist.fct: sumofsquares")

    expect_identical(X[["centers"]], 10L)
    X[["centers"]] <- 2L
    expect_identical(X[["centers"]], 2L)

    X <- SOMParam(centers=log)
    expect_true(is.function(X[["centers"]]))

    X[["radius"]] <- c(0, 0.1)
    expect_output(show(X), "radius: 0 0.1")

    X <- SOMParam(centers=10, rlen=50)
    expect_identical(X[["rlen"]], 50L)
    X[["rlen"]] <- 100L 
    expect_identical(X[["rlen"]], 100L)
})

test_that("clusterRows works correctly", {
    m <- matrix(runif(10000), ncol=10)

    set.seed(100)
    out <- clusterRows(m, SOMParam(25))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 25L)

    # Beahves the same with a seed.
    set.seed(100)
    out2 <- clusterRows(m, SOMParam(25))
    expect_identical(out, out2)

    set.seed(200)
    out3 <- clusterRows(m, SOMParam(25))
    expect_false(identical(out, out3))

    # Consistent with reference values.
    set.seed(100)
    ref <- kohonen::som(m, kohonen::somgrid(5,5))
    expect_identical(out, factor(ref$unit.classif))

    # Trying with fewer cells.
    m <- matrix(runif(1000), ncol=10)
    colnames(m) <- seq_len(ncol(m))

    set.seed(100)
    out <- clusterRows(m, SOMParam(12, dim.ratio=3/4))
    expect_true(is.factor(out))
    expect_identical(length(out), nrow(m))
    expect_identical(nlevels(out), 12L)

    set.seed(100)
    ref <- kohonen::som(m, kohonen::somgrid(3,4))
    expect_identical(out, factor(ref$unit.classif))
})

test_that("clusterRows responds to the functions and full=TRUE", {
    m <- matrix(runif(12960), ncol=10)

    set.seed(9999)
    out <- clusterRows(m, SOMParam(sqrt))
    expect_identical(length(out), nrow(m))
    expect_equal(nlevels(out), round(sqrt(nrow(m))))

    set.seed(9999)
    full <- clusterRows(m, SOMParam(sqrt), full=TRUE)
    expect_identical(out, full$cluster)
    expect_true(is.list(full$objects))
    expect_identical(names(full$objects), "som")
})

test_that("clusterRows responds to the options", {
    m <- matrix(runif(10000), ncol=10)

    set.seed(100000)
    suppressWarnings(ref<- clusterRows(m, SOMParam(16)))

    set.seed(100000)
    suppressWarnings(con <- clusterRows(m, SOMParam(16)))
    expect_true(identical(ref, con)) # as a control

    set.seed(100000)
    suppressWarnings(out <- clusterRows(m, SOMParam(16, rlen=20, alpha=c(0.12, 0.06)))) 
    expect_false(identical(ref, out)) 

    set.seed(100000)
    ref <- kohonen::som(m, kohonen::somgrid(4,4), rlen=20, alpha=c(0.12, 0.06))
    expect_identical(out, factor(ref$unit.classif))
})
