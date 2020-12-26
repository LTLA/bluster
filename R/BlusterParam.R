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

