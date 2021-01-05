setValidity("FixedNumberParam", function(object) {
    msg <- character(0)

    centerx <- centers(object)
    if (is.integer(centerx)) {
        if (!.positive_number(centerx)) {
            msg <- c(msg, "integer 'centers' must be a positive scalar")
        }
    }

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

#' @export
setMethod("centers", "FixedNumberParam", function(x, n=NULL) {
    xcenters <- x@centers
    if (is.integer(xcenters) || is.null(n)) {
        xcenters
    } else {
        as.integer(round(xcenters(n)))
    }
})

#' @export
setReplaceMethod("centers", "FixedNumberParam", function(x, value) {
    x@centers <- value
    x
})

#' @export
setMethod("show", "FixedNumberParam", function(object) {
    callNextMethod()
    cat(sprintf("centers: %s\n", if (is.function(object@centers)) "variable" else object@centers))
})
