\name{LogrankTest}
\alias{logrank_test}
\alias{logrank_test.formula}
\alias{logrank_test.IndependenceProblem}
\title{ Independent Two- and K-Sample Tests for Censored Data }
\description{
    
}
\usage{
\method{logrank_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{logrank_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...)
}
\arguments{
  \item{formula}{a formula of the form \code{Surv(time, event) ~ x | block} where 
    \code{time} is a positive numeric variable denoting the survival time and
    \code{event} is a logical being \code{TRUE} when the event of interest
    was observed and \code{FALSE} in case of censoring. \code{x} is a factor
    with two or more levels giving the corresponding groups. \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
      integer valued weights for the observations.}
  \item{x}{an object of class \code{IndependenceProblem}.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of the survival
  functions in the groups induced by \code{x} is tested. 

}
\value{
  An object inheriting from class \code{ScalarIndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}. 

}
\examples{

data(ocarcinoma)

logrank_test(Surv(time, cens) ~ stadium, data = ocarcinoma, 
             distribution = "exact")

}
\keyword{htest}
