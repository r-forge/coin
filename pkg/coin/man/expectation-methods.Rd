\name{expectation-methods}
\docType{methods}
\alias{expectation}
\alias{expectation-methods}
\alias{expectation,IndependenceLinearStatistic-method}
\alias{expectation,IndependenceTest-method}
\alias{variance}
\alias{variance-methods}
\alias{variance,Variance-method}
\alias{variance,CovarianceMatrix-method}
\alias{variance,IndependenceLinearStatistic-method}
\alias{variance,IndependenceTest-method}
\alias{covariance}
\alias{covariance-methods}
\alias{covariance,CovarianceMatrix-method}
\alias{covariance,IndependenceLinearStatistic-method}
\alias{covariance,QuadTypeIndependenceTestStatistic-method}
\alias{covariance,IndependenceTest-method}
\title{Extraction of the Expectation, Variance and Covariance of the Linear
  Statistic}
\description{
  Methods for extraction of the expectation, variance and covariance of the
  linear statistic.
}
\usage{
\S4method{expectation}{IndependenceLinearStatistic}(object, \dots)
\S4method{expectation}{IndependenceTest}(object, \dots)

\S4method{variance}{Variance}(object, \dots)
\S4method{variance}{CovarianceMatrix}(object, \dots)
\S4method{variance}{IndependenceLinearStatistic}(object, \dots)
\S4method{variance}{IndependenceTest}(object, \dots)

\S4method{covariance}{CovarianceMatrix}(object, \dots)
\S4method{covariance}{IndependenceLinearStatistic}(object, invert = FALSE, \dots)
\S4method{covariance}{QuadTypeIndependenceTestStatistic}(object, invert = FALSE, \dots)
\S4method{covariance}{IndependenceTest}(object, invert = FALSE, \dots)
}
\arguments{
  \item{object}{
    an object from which the expectation, variance or covariance of the linear
    statistic can be extracted.
  }
  \item{invert}{
    a logical indicating that the Moore-Penrose inverse of the covariance should
    be computed.  Defaults to \code{FALSE}.
  }
  \item{\dots}{
    further arguments (currently ignored).
  }
}
\details{
  The methods \code{expectation}, \code{variance} and \code{covariance} extract
  the expectation, variance and covariance, respectively, of the linear
  statistic.
}
\value{
  The expectation, variance or covariance of the linear statistic extracted from
  \code{object}.  A matrix.
}
\examples{\dontshow{suppressWarnings(RNGversion("3.5.2")); set.seed(711109)}
## Example data
dta <- data.frame(
    y = gl(3, 2),
    x = sample(gl(3, 2))
)

## Asymptotic Cochran-Mantel-Haenszel Test
ct <- cmh_test(y ~ x, data = dta)

## The linear statistic, i.e., the contingency table...
(T <- statistic(ct, type = "linear"))

## ...and its expectation...
(mu <- expectation(ct))

## ...and covariance...
(Sigma <- covariance(ct))

## ...and its inverse
(SigmaPlus <- covariance(ct, invert = TRUE))

## The standardized contingency table...
(T - mu) / sqrt(variance(ct))

## ...is identical to the standardized linear statistic
statistic(ct, type = "standardized")

## The quadratic form...
U <- as.vector(T - mu)
U \%*\% SigmaPlus \%*\% U

## ...is identical to the test statistic
statistic(ct, type = "test")
}
\keyword{methods}
