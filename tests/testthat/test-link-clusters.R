# This tests the linkClusters function.
# library(testthat); library(bluster); source("test-link-clusters.R")

set.seed(9999999)
x1 <- sample(LETTERS[1:10], 100, replace=TRUE)
x2 <- sample(1:15, 100, replace=TRUE)
x3 <- sample(letters[1:5], 100, replace=TRUE)
together <- list(first=x1, second=x2, third=x3)

set.seed(99999999)
test_that("linkClusters works correctly", {
    linked <- linkClusters(together, denominator="min")
    ulinked <- linkClusters(together, denominator="union")
    xlinked <- linkClusters(together, denominator="max")

    # Prefixing works correctly.
    out <- names(igraph::V(linked))
    expect_true(all(grepl("^(first|second|third)", out)))

    edges <- igraph::ends(linked, igraph::E(linked))
    expect_identical(nrow(edges), length(unique(c(paste(x1, x2), paste(x1, x3), paste(x2, x3)))))

    expect_identical(edges, igraph::ends(ulinked, igraph::E(ulinked)))
    expect_identical(edges, igraph::ends(xlinked, igraph::E(xlinked)))

    # Edge weight calculations are performed correctly.
    ref.weights <- uref.weights <- xref.weights <- numeric(nrow(edges)) 
    for (e in seq_len(nrow(edges))) {
        left <- edges[e,1]
        left.cluster <- sub("\\..*", "", left)
        left.level <- sub(".*\\.", "", left)

        right <- edges[e,2]
        right.cluster <- sub("\\..*", "", right)
        right.level <- sub(".*\\.", "", right)

        has.left <- together[[left.cluster]]==left.level 
        has.right <- together[[right.cluster]]==right.level 
        shared <- sum(has.left & has.right)
        n.left <- sum(has.left)
        n.right <- sum(has.right)

        ref.weights[e] <- shared / min(n.left, n.right)
        uref.weights[e] <- shared / sum(has.left | has.right)
        xref.weights[e] <- shared / max(n.left, n.right)
    }

    expect_identical(ref.weights, igraph::E(linked)$weight)
    expect_identical(uref.weights, igraph::E(ulinked)$weight)
    expect_identical(xref.weights, igraph::E(xlinked)$weight)
})

set.seed(999999991)
test_that("linkClusters behaves correctly with unnumbered elements", {
    ref <- linkClusters(together)
    stripped <- ref[]
    nofix <- sub(".*\\.", "", rownames(stripped))
    colnames(stripped) <- nofix
    rownames(stripped) <- nofix

    linked <- linkClusters(together, prefix=FALSE)
    expect_identical(stripped, linked[])

    ref <- linkClusters(unname(together))
    together2 <- together
    names(together2) <- seq_along(together)
    linked <- linkClusters(together2)
    expect_identical(linked[], ref[])
})

test_that("linkClusters handles error states smoothly", {
    expect_error(linkClusters(list(letters, 1:10)), "same length")

    empty <- linkClusters(list(integer(0), character(0)))
    expect_identical(length(igraph::E(empty)), 0L)
    expect_identical(length(igraph::V(empty)), 0L)
})
