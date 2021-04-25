#' @export
setMethod("[[", "BlusterParam", function(x, i) {
    slot(x, i)
})

#' @export
setReplaceMethod("[[", "BlusterParam", function(x, i, j, ..., value) {
    slot(x, i) <- value
    x
})

#' @export
setMethod("show", "BlusterParam", function(object) {
    cat(paste0("class: ", class(object)[1], "\n"))
})

setValidity2("BlusterParam", function(object) {
    what <- .defaultScalarArguments(object)
    for (x in names(what)) {
        val <- slot(object, x)
        if (!is.null(val)) {
            if (length(val)!=1 || !is(val, what[[x]])) {
                return(sprintf("'%s' should be NULL or a %s scalar", x, what[[x]]))
            }
        }
    }
})

#' @export
setMethod(".defaultScalarArguments", "BlusterParam", function(x) character(0))

#' @export
#' @rdname defaultArguments
.showScalarArguments <- function(object) {
    what <- .defaultScalarArguments(object)
    for (x in names(what)) {
        val <- slot(object, x)
        if (is.null(val)) {
            val <- "[default]"
        }
        cat(sprintf("%s: %s\n", x, val))
    }
}

#' @export
#' @rdname defaultArguments
.extractScalarArguments <- function(x) {
    args <- list()
    what <- .defaultScalarArguments(x)
    for (i in names(what)) {
        val <- slot(x, i)
        if (!is.null(val)) {
            args[[i]] <- val
        }
    }
    args

}
