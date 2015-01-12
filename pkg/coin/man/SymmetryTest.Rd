\name{SymmetryTest}
\alias{symmetry_test}
\alias{symmetry_test.formula}
\alias{symmetry_test.table}
\alias{symmetry_test.SymmetryProblem}

\title{ General Symmetry Test }

\description{
  Testing the symmetry of a response for repeated measurements in a
  complete block design.
}
\usage{
\method{symmetry_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{symmetry_test}{table}(object, \dots)
\method{symmetry_test}{SymmetryProblem}(object,
    teststat = c("maximum", "quadratic", "scalar"),
    distribution = c("asymptotic", "approximate", "exact", "none"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL,
    check = NULL, paired = FALSE, \dots)
}
\arguments{
  \item{formula}{a formula of the form \code{y1 + ... + yq ~ x | block} where
  \code{y1} to \code{yq} may be measured on arbitrary scales, \code{x} is a
  factor with two or more levels giving the corresponding groups and
  \code{block} is an optional factor (which is generated automatically when
  omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}
  \item{object}{an object inheriting from class \code{"SymmetryProblem"}
    or an object of class \code{"table"}.}
  \item{teststat}{a character, the type of test statistic to be applied:
    either a standardized scalar test statistic (\code{"scalar"}), or a
    maximum type statistic (\code{"maximum"}, default) or a quadratic form
    (\code{"quadratic"}).}
  \item{distribution}{a character, the null distribution of the test
    statistic can be computed \code{"exact"}ly or can be approximated by
    its asymptotic distribution (\code{"asymptotic"}, default) or via
    Monte Carlo resampling (\code{"approximate"}).  Alternatively, the
    functions \code{\link{exact}}, \code{\link{approximate}} or
    \code{\link{asymptotic}} can be used to specify how the exact
    conditional distribution of the test statistic should be calculated
    or approximated.  It is also possible to specify a function with one
    argument (taking objects inheriting from
    \code{"IndependenceTestStatistic"}) which return an object of class
    \code{"NullDistribution"}.  Computation of the null distribution can
    be suppressed by specifying \code{"none"}.}
  \item{alternative}{a character, the alternative hypothesis must be one
    of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.}
  \item{xtrafo}{a function of transformations (see \code{\link{trafo}})
    to be applied to the variables on the right hand side of
    \code{formula}, see below.}
  \item{ytrafo}{a function of transformations (see \code{\link{trafo}})
    to be applied to the variables on the left hand side of
    \code{formula}, see below.}
  \item{scores}{a named list of scores to be attached to ordered
    factors.  In case a variable is an unordered factor, it is coerced
    to \code{ordered} first.}
  \item{check}{a function to be applied to objects of class
    \code{"SymmetryTest"} in order to check for specific forms, see
    below.}
  \item{paired}{a logical, indicating that paired data have been
    transformed in such a way that the resulting unstandardized linear
    statistic is the sum of the absolute values of the positive
    differences between the paired observations.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  This function handles symmetry problems where each block corresponds to
  exactly one observation with repeated measurements.

}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\examples{
## permutation test for paired observations
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
dta <- data.frame(
    y = c(y1, y2),
    x = gl(2, length(y1)),
    block = factor(rep(seq_along(y1), 2))
)
symmetry_test(y ~ x | block, data = dta,
              distribution = "exact", alternative = "greater")

## alternatively: transform data and set 'paired = TRUE'
diff <- y1 - y2
y <- as.vector(rbind(abs(diff) * (diff >= 0), abs(diff) * (diff < 0)))
x <- factor(rep(0:1, length(diff)), labels = c("pos", "neg"))
block <- gl(length(diff), 2)
symmetry_test(y ~ x | block,
              distribution = "exact", alternative = "greater",
              paired = TRUE)
}
\keyword{htest}
