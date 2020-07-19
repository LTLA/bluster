.non_na_scalar <- function(val) {
    length(val)==1L && !is.na(val)
}

.non_na_number <- function(val) {
    .non_na_scalar(val) && is.numeric(val)
}

.positive_number <- function(val) {
    .non_na_number(val) && val > 0
}
