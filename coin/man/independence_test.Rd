\name{IndependenceTest}
\alias{independence_test}
\title{ General Independence Test } 
\description{
    The independence between two sets of variables of arbitrary measurement
    scales, possibly stratified in blocks, is tested conditional on the data.
}
\usage{
\method{independence_test}{formula}(formula, data, subset = NULL, \dots)
\method{independence_test}{IndependenceProblem}(x, 
    teststat = c("maxtype", "quadtype", "scalar"),
    distribution = c("asympt", "approx", "exact"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, check = NULL, ...)
\method{independence_test}{table}(x, 
    distribution = c("asympt", "approx"), ...)
}
\arguments{
  \item{formula}{a formula of the form 
    \code{y1 + y2 + ... + yp ~ x1 + x2 + ... xq | block} where the variables
    on the left and right hand side may be measured on arbitrary scales
    (including censored ones) and \code{block} is an optional factor for 
    stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{x}{an object inheriting from class \code{IndependenceProblem} or an
           object of class \code{table}.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{xtrafo}{a function of transformations to be applied to the right hand side of
    \code{formula} or \code{x@x} see below.}
  \item{ytrafo}{a function of transformations to be applied to the left hand side of
    \code{formula} or \code{x@y} see below.}
  \item{check}{a function to be applied to objects of class
      \code{IndendenceTest} in order to check for specific forms, see below.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the independence between the left hand side and the
  right hand side of \code{formula}, possibly stratified by \code{block}, is
  tested. 

}
\value{
  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.
}
\examples{

data(bloodp)

### independence of bp and group
independence_test(bp ~ group, data = bloodp,

    ### exact null distribution
    distribution = "exact", 

    ### one-sided test
    alternative = "less",

    ### apply normal scores to bp
    ytrafo = function(data) trafo(data, numeric_trafo = normal_trafo),

    ### indicator matrix of 2nd level of group
    xtrafo = function(data) trafo(data, factor_trafo = function(x)
        matrix(x == levels(x)[2], ncol = 1)),

    ### check if this is a independent two sample problem
    check = is_2sample)

### same as
normal_test(bp ~ group, data = bloodp, distribution = "exact", 
            alternative = "greater")
}
\keyword{htest}
