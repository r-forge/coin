\name{SpearmanTest}
\alias{spearman_test.formula}
\alias{spearman_test.IndependenceProblem}
\alias{spearman_test}
\title{ Spearman's Test on Independence }
\description{
    Testing the independence of two numeric variables.
}
\usage{
\method{spearman_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{spearman_test}{IndependenceProblem}(x,
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx"), ...) 
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    and \code{x} are numeric variables and \code{block} is an
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
    can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the independence of \code{y} and \code{x} 
  is tested. 

}
\value{
  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\examples{

spearman_test(CONT ~ INTG, data = USJudgeRatings)

}
\keyword{htest}
