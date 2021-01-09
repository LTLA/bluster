.non_na_scalar <- function(val) {
    length(val)==1L && !is.na(val)
}

.non_na_number <- function(val) {
    .non_na_scalar(val) && is.numeric(val)
}

.positive_number <- function(val) {
    .non_na_number(val) && val > 0
}

.check_positive_slots <- function(object, names) {
    msg <- character(0)
    for (i in names) {
        val <- slot(object, i)
        if (!is.null(val) && !.positive_number(val)) {
            msg <- c(msg, sprintf("'%s' should be positive", i))
        }
    }
    msg
}

.check_nonna_slots <- function(object, names) {
    msg <- character(0)
    for (i in names) {
        val <- slot(object, i)
        if (!is.null(val) && !.non_na_scalar(val)) {
            msg <- c(msg, sprintf("'%s' should not be NA", i))
        }
    }
    msg
}
