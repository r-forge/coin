\name{ExpectationCovarianceStatistic}
\alias{ExpectationCovarianceStatistic}
\title{
  \pkg{libcoin} C API Example
}
\description{
  A working example of calling compiled code from the \pkg{libcoin} package.
}
\usage{
ExpectationCovarianceStatistic(x, y)
}
\arguments{
  \item{x}{A numeric matrix.}
  \item{y}{A numeric matrix.}
}
\details{
  This function calls the C function  \code{R_ExpectationCovarianceStatistic}, a
  simple wrapper for the C function
  \code{libcoin_R_ExpectationCovarianceStatistic} defined in the \pkg{libcoin}
  package.

  For more information on this approach, consult section 5.4.2 of the Writing R
  Extensions manual.
}
\value{
  A list.
}
\examples{
n <- 100
p <- 4
q <- 2
X <- matrix(runif(p * n), nc = p)
Y <- matrix(runif(q * n), nc = q)
ExpectationCovarianceStatistic(X, Y)
}
