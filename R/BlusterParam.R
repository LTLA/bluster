#' @export
setMethod("[[", "BlusterParam", function(x, i) {
    if (i %in% slotNames(x)) {
        slot(x, i)
    } else {
        slot(x, .extras(x))[[i]]
    }
})

#' @export
setReplaceMethod("[[", "BlusterParam", function(x, i, j, ..., value) {
    if (i %in% slotNames(x)) {
        slot(x, i) <- value
    } else {
        ex.name <- .extras(x)
        extras <- slot(x, ex.name)
        extras[[i]] <- value
        slot(x, ex.name) <- extras
    }
    x
})

#' @export
setMethod("show", "BlusterParam", function(object) {
    cat(paste0("class: ", class(object)[1], "\n"))
})

