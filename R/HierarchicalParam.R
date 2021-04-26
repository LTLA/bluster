#' @export
setMethod(".defaultScalarArguments", "HierarchicalParam", function(x) c(metric="character"))

#' @importFrom S4Vectors setValidity2
setValidity2("HierarchicalParam", function(object) {
    if (!.non_na_scalar(slot(object, "cut.dynamic"))) {
        return("'cut.dynamic' must be a non-missing logical scalar")
    }
    TRUE
})

#' @importFrom S4Vectors coolcat
setMethod("show", "HierarchicalParam", function(object) {
    callNextMethod()
    .showScalarArguments(object)

    fun <- object@cut.fun
    if (is.null(fun)) {
        if (object@cut.dynamic) {
            cat("cut.fun: cutreeDynamic\n")
        } else {
            cat("cut.fun: cutree\n")
        }
    } else {
        cat("cut.fun: custom\n")
    }
    coolcat("cut.params(%i): %s", names(object@cut.params))
})

.cut_hierarchical <- function(hcl, dst, BLUSPARAM) {
    fun <- BLUSPARAM@cut.fun
    args <- BLUSPARAM@cut.params

    if (is.null(fun)) {
        if (!BLUSPARAM@cut.dynamic) {
            fun <- cutree
            if (is.null(args$k)) {
                if (is.null(args$h)) {
                    args$h <- max(hcl$height)/2
                }
            }

        } else {
            fun <- function(...) unname(dynamicTreeCut::cutreeDynamic(..., verbose=0))
            args$dist <- as.matrix(dst)
        }
    }

    clusters <- do.call(fun, c(list(hcl), args))
    factor(clusters)
}
