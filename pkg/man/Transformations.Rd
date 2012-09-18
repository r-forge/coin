\name{Transformations}
\alias{trafo}
\alias{id_trafo}
\alias{rank_trafo}
\alias{ansari_trafo}
\alias{fligner_trafo}
\alias{normal_trafo}
\alias{median_trafo}
\alias{consal_trafo}
\alias{maxstat_trafo}
\alias{fmaxstat_trafo}
\alias{logrank_trafo}
\alias{f_trafo}
\alias{of_trafo}
\alias{mcp_trafo}
\title{ Functions for Data Transformations }
\description{
  Rank-transformations for numerical data or dummy codings of factors.
}
\usage{
trafo(data, numeric_trafo = id_trafo, factor_trafo = f_trafo,
    ordered_trafo = of_trafo, surv_trafo = logrank_trafo,
    var_trafo = NULL, block = NULL)
id_trafo(x)
rank_trafo(x, ties.method = c("mid-ranks", "random"))
ansari_trafo(x, ties.method = c("mid-ranks", "average-scores"))
fligner_trafo(x, ties.method = c("mid-ranks", "average-scores"))
normal_trafo(x, ties.method = c("mid-ranks", "average-scores"))
median_trafo(x, mid.score = c("0", "0.5", "1"))
consal_trafo(x, ties.method = c("mid-ranks", "average-scores"), a = 5)
maxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
fmaxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
logrank_trafo(x, ties.method = c("logrank", "HL", "average-scores"))
f_trafo(x)
of_trafo(x)
mcp_trafo(...)
}
\arguments{
  \item{data}{an object of class \code{data.frame}.}
  \item{numeric_trafo}{a function to by applied to \code{numeric}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns.}
  \item{factor_trafo}{a function to by applied to \code{factor}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually a dummy or contrast
    matrix).}
  \item{ordered_trafo}{a function to by applied to \code{ordered}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually some scores).}
  \item{surv_trafo}{a function to by applied to
    elements of class \code{Surv} of \code{data} returning a
    matrix with \code{nrow(data)} rows and an arbitrary number of columns.}
  \item{var_trafo}{an optional named list of functions to be applied to the
                   corresponding variables in \code{data}.}
  \item{block}{an optional factor those levels are interpreted as blocks.
    \code{trafo} is applied to each level of \code{block} separately.}
  \item{x}{an object of classes \code{numeric}, \code{ordered}, \code{factor} or
    \code{Surv}.}
  \item{ties.method}{method used to adjust the scores for ties. Either
    the score generating function is applied to \code{"mid-ranks"} or
    scores, based on random ranks, are averaged
    (\code{"average-scores"}). When using \code{rank_trafo}, ties can be
    broken at \code{"random"}. For ties handling in case of censored
    data, see \code{\link{surv_test}}.}
  \item{mid.score}{a character, the score assigned to observations
    exactly equal to the median.  Must be one of \code{"0"} (default),
    \code{"0.5"} or \code{"1"}; \code{"0.5"} leads to symmetric median
    scores ensuring that the resulting test statistic has a symmetric
    distribution, or at least approximately so, under the null
    hypothesis.}
  \item{a}{a numeric vector, the value taken as the constant \eqn{a}
    in the Conover-Salsburg scores.  Defaults to 5.}
  \item{minprob}{a fraction between 0 and 0.5.}
  \item{maxprob}{a fraction between 0.5 and 1.}
  \item{\dots}{factor name and contrast matrix (as matrix or character)
               in a tag = value format for multiple comparisons
               for one single unordered factor, see \code{\link[multcomp]{mcp}}.}
}
\details{

  The utility functions documented here are used to define special
  independence tests.

  \code{trafo} applies its arguments to the elements of \code{data}
  according to the classes of the elements.

  \code{id_trafo} is the identity transformation and \code{f_trafo}
  computes dummy matrices for factors.

  \code{rank_trafo} is the rank transformation.

  \code{ansari_trafo} and \code{fligner_trafo} compute Ansari-Bradley
  or Fligner scores for scale problems.

  \code{normal_trafo}, \code{median_trafo} and \code{consal_trafo}
  implement normal scores, median scores or Conover-Salsburg scores (see
  \code{\link{neuropathy}}) for location problems, \code{logrank_trafo}
  returns logrank scores for censored data.

  A \code{trafo} function with modified default arguments is usually
  feeded into \code{\link{independence_test}} via the \code{xtrafo}
  or \code{ytrafo} arguments.

  Fine tuning (different transformations for different variables) is
  possible by supplying a named list of functions to the \code{var_trafo}
  argument.

}
\note{

  In \pkg{coin} version 1.1-0 (and later), all transformation functions
  are now passing through missing values (i.e., \code{NA}s).
  Furthermore, \code{median_trafo} and \code{logrank_trafo} are now
  increasing functions of the data (in conformity with most other
  transformations in \pkg{coin}).

}
\value{

  A named matrix with \code{nrow(data)} rows and
  arbitrary number of columns. User-supplied transformations must
  return a numeric vector or matrix.

}
\examples{

  ### dummy matrices, 2-sample problem (only one column)
  f_trafo(y <- gl(2, 5))

  ### score matrices
  of_trafo(y <- ordered(gl(3, 5)))

  ### K-sample problem (K columns)
  f_trafo(y <- gl(5, 2))

  ### normal scores
  normal_trafo(x <- rnorm(10))

  ### and now together
  trafo(data.frame(x = x, y = y), numeric_trafo = normal_trafo)

  ### the same, more flexible when multiple variables are in play
  trafo(data.frame(x = x, y = y), var_trafo = list(x = normal_trafo))

  ### maximally selected statistics
  maxstat_trafo(rnorm(10))

  ### apply transformation blockwise (e.g. for Friedman test)
  trafo(data.frame(y = 1:20), numeric_trafo = rank_trafo, block = gl(4, 5))

  ### multiple comparisons
  mcp_trafo(x = "Tukey")(data.frame(x = gl(3, 3)))

}
\keyword{manip}
