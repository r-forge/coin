
ExpectationCovarianceStatistic <- function(x, y) {
    .Call(R_ExpectationCovarianceStatistic,
          x, y, integer(0), integer(0), integer(0), 0L,
          sqrt(.Machine$double.eps))
}
