\name{MaxstatTest}
\alias{maxstat_test.formula}
\alias{maxstat_test.IndependenceProblem}
\alias{maxstat_test}
\title{ Maximally selected test statistics }
\description{
    Testing the independence of a set of ordered or numeric covariates and a
    response of arbitrary measurement scale against cutpoint alternatives. 
}
\usage{
\method{maxstat_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{maxstat_test}{IndependenceProblem}(x, 
    distribution = c("asympt", "approx", ...)
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x1 + ... + xp | block} where \code{y}
    is a variable measured at arbitrary scale and the covariates \code{x1}
    to \code{xp} are at least of class \code{ordered}; \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}   
  \item{x}{an object inheriting from class \code{IndependenceProblem}.}
  \item{distribution}{the null distribution of the test statistic
    can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of independence of all covariates to the response
  \code{y} against cutpoint alternatives is tested.

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\examples{

data(treepipit)

maxstat_test(counts ~ coverstorey, data = treepipit)

}
\keyword{htest}
