\name{Distribution}
\alias{exact}
\alias{asymptotic}
\alias{approximate}

\title{ Distribution under the Null Hypothesis }

\description{
  Specification of the exact, or approximation of the exact, conditional
  distribution of test statistics under the null hypothesis.
}
\usage{
exact(algorithm = c("auto", "shift", "split-up"), fact = NULL)
approximate(B = 10000)
asymptotic(maxpts = 25000, abseps = 0.001, releps = 0)
}
\arguments{
  \item{algorithm}{a character, specifying the algorithm to be used for the
  computation of the exact conditional distribution: either \code{"auto"}
  (default), \code{"shift"} or \code{"split-up"}.}
  \item{fact}{a positive integer to multiply the response values with.}
  \item{B}{a positive integer, the number of Monte Carlo replications to
    approximate the exact conditional distribution.}
  \item{maxpts}{a positive integer, the maximum number of function
    values, see \code{\link[mvtnorm]{pmvnorm}}.}
  \item{abseps}{double, the absolute error tolerance, see
    \code{\link[mvtnorm]{pmvnorm}}.}
  \item{releps}{double, the relative error tolerance, see
    \code{\link[mvtnorm]{pmvnorm}}.}
}
\details{

  The \code{distribution} argument to \code{\link{independence_test}}
  can be specified with additional arguments using those functions.

  Exact algorithms are currently only implemented for two-sample problems.

}
\keyword{htest}
