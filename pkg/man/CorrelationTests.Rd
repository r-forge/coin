\name{CorrelationTests}
\alias{spearman_test}
\alias{spearman_test.formula}
\alias{spearman_test.IndependenceProblem}
\alias{fisyat_test}
\alias{fisyat_test.formula}
\alias{fisyat_test.IndependenceProblem}
\alias{quadrant_test}
\alias{quadrant_test.formula}
\alias{quadrant_test.IndependenceProblem}
\alias{koziol_test}
\alias{koziol_test.formula}
\alias{koziol_test.IndependenceProblem}

\encoding{UTF-8}

\title{ Correlation Tests }

\description{
  Testing the independence of two numeric variables.
}
\usage{

\method{spearman_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{spearman_test}{IndependenceProblem}(object,
    distribution = c("asymptotic", "approximate"), \dots)

\method{fisyat_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{fisyat_test}{IndependenceProblem}(object,
    distribution = c("asymptotic", "approximate"),
    ties.method = c("mid-ranks", "average-scores"), \dots)

\method{quadrant_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{quadrant_test}{IndependenceProblem}(object,
    distribution = c("asymptotic", "approximate"),
    mid.score = c("0", "0.5", "1"), \dots)

\method{koziol_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{koziol_test}{IndependenceProblem}(object,
    distribution = c("asymptotic", "approximate"),
    ties.method = c("mid-ranks", "average-scores"), \dots)

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
  \item{object}{an object of class \code{IndependenceProblem}.}
  \item{distribution}{a character, the null distribution of the test
    statistic can be approximated by its asymptotic distribution
    (\code{"asymptotic"}, default) or via Monte Carlo resampling
    (\code{"approximate"}).  Alternatively, the functions
    \code{\link{approximate}} or \code{\link{asymptotic}} can be used to
    specify how the exact conditional distribution of the test statistic
    should be approximated.}
  \item{ties.method}{a character, two methods are available to adjust
    scores for ties, either the score generating function is applied to
    \code{"mid-ranks"} (default) or the scores computed based on random
    ranks are averaged for all tied values (\code{"average-scores"}).}
  \item{mid.score}{a character, the score assigned to observations
    exactly equal to the median.  Must be one of \code{"0"} (default),
    \code{"0.5"} or \code{"1"}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the independence of \code{y} and \code{x}
  is tested.

  The test procedures documented here perform the Spearman correlation
  test (\code{spearman_test}), the Fisher-Yates correlation test using
  van der Waerden scores (\code{fisyat_test}), the quadrant test
  (\code{quadrant_test}), and the Koziol-Nemec test
  (\code{koziol_test}).  For a general description of these tests we
  refer to Hajek, Sidak and Sen (1999, Section 4.6).  The Koziol-Nemec
  test was suggested by Koziol & Nemec (1979).

  For the adjustment of scores for tied values see Hajek, Sidak and Sen
  (1999, p. 133--135).

}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
 }
\references{

  Jaroslav Hajek, Zbynek Sidak & Pranab K. Sen (1999).
  \emph{Theory of Rank Tests}, 2nd Edition.
  San Diego: Academic Press.

  James A. Koziol & Amanda F. Nemec (1979).
  On a \enc{Cramér}{Cramer}-von Mises type statistic for testing
  bivariate independence.
  \emph{The Canadian Journal of Statistics} \bold{7}(1), 43--52.

}
\examples{

  spearman_test(CONT ~ INTG, data = USJudgeRatings)

  fisyat_test(CONT ~ INTG, data = USJudgeRatings)

  quadrant_test(CONT ~ INTG, data = USJudgeRatings)

  koziol_test(CONT ~ INTG, data = USJudgeRatings)

}
\keyword{htest}
