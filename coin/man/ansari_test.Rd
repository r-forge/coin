\name{ScaleTests}
\alias{ansari_test}
\alias{fligner_test}
\title{ Independent Two- and K-Sample Scale Tests }
\description{
    
}
\usage{
\method{ansari_test}{formula}(formula, data, subset = NULL, \dots)
\method{ansari_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)
\method{fligner_test}{formula}(formula, data, subset = NULL, \dots)
\method{fligner_test}{IndependenceProblem}(x, 
    distribution = c("asympt", "approx"), ...)
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric variable giving the data values and \code{x} a factor
    with two or more levels giving the corresponding groups. \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{x}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{conf.int}{a logical indicating whether a confidence interval
    for the difference in location should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested. In particular, the methods
  documented here are designed to detect scale alternatives. 

}
\value{
  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue}, \code{statistic} and (optionally)
  \code{confint}.

}
\examples{

data(bloodp)

# Ansari-Bradley test

ansari_test(bp ~ group, data = bloodp)
ansari_test(bp ~ group, data = bloodp, distribution = "exact")
ansari_test(bp ~ group, data = bloodp, distribution = "approx", B = 9999)

}
\keyword{htest}
