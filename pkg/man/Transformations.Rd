\name{Transformations}
\alias{trafo}
\alias{id_trafo}
\alias{rank_trafo}
\alias{normal_trafo}
\alias{median_trafo}
\alias{savage_trafo}
\alias{consal_trafo}
\alias{koziol_trafo}
\alias{ansari_trafo}
\alias{fligner_trafo}
\alias{logrank_trafo}
\alias{logrank_weight}
\alias{maxstat_trafo}
\alias{fmaxstat_trafo}
\alias{ofmaxstat_trafo}
\alias{f_trafo}
\alias{of_trafo}
\alias{mcp_trafo}

\title{ Functions for Data Transformations }

\description{
  Transformations for numerical data or dummy codings of factors.
}
\usage{
trafo(data, numeric_trafo = id_trafo, factor_trafo = f_trafo,
    ordered_trafo = of_trafo, surv_trafo = logrank_trafo,
    var_trafo = NULL, block = NULL)
id_trafo(x)
rank_trafo(x, ties.method = c("mid-ranks", "random"))
normal_trafo(x, ties.method = c("mid-ranks", "average-scores"))
median_trafo(x, mid.score = c("0", "0.5", "1"))
savage_trafo(x, ties.method = c("mid-ranks", "average-scores"))
consal_trafo(x, ties.method = c("mid-ranks", "average-scores"), a = 5)
koziol_trafo(x, ties.method = c("mid-ranks", "average-scores"), j = 1)
ansari_trafo(x, ties.method = c("mid-ranks", "average-scores"))
fligner_trafo(x, ties.method = c("mid-ranks", "average-scores"))
logrank_trafo(x, ties.method = c("mid-ranks", "Hothorn-Lausen", "average-scores"),
              weight = logrank_weight, \dots)
logrank_weight(time, n.risk, n.event,
               type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Prentice",
                        "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
                        "Fleming-Harrington", "Self"),
               rho = NULL, gamma = NULL)
maxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
fmaxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
ofmaxstat_trafo(x, minprob = 0.1, maxprob = 1 - minprob)
f_trafo(x)
of_trafo(x, scores = NULL)
mcp_trafo(\dots)
}
\arguments{
  \item{data}{an object of class \code{data.frame}.}
  \item{numeric_trafo}{a function to be applied to \code{numeric}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns.}
  \item{factor_trafo}{a function to be applied to \code{factor} elements
    of \code{data} returning a matrix with \code{nrow(data)} rows and an
    arbitrary number of columns (usually a dummy or contrast matrix).}
  \item{ordered_trafo}{a function to be applied to \code{ordered}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually some scores).}
  \item{surv_trafo}{a function to be applied to elements of class
    \code{Surv} of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns.}
  \item{var_trafo}{an optional named list of functions to be applied to
    the corresponding variables in \code{data}.}
  \item{block}{an optional factor those levels are interpreted as
    blocks.  \code{trafo} is applied to each level of \code{block}
    separately.}
  \item{x}{an object of classes \code{numeric}, \code{ordered},
    \code{factor} or \code{Surv}.}
  \item{ties.method}{method used to adjust the scores for ties.  Either
    the score generating function is applied to \code{"mid-ranks"} or
    scores, based on random ranks, are averaged
    (\code{"average-scores"}).  When using \code{rank_trafo}, ties can
    be broken at \code{"random"}.  For ties handling in the case of
    censored data, see \code{\link{surv_test}}.}
  \item{mid.score}{a character, the score assigned to observations
    exactly equal to the median.  Must be one of \code{"0"} (default),
    \code{"0.5"} or \code{"1"}; \code{"0.5"} leads to symmetric median
    scores (see also \code{\link{median_test}}).}
  \item{a}{a numeric vector, the value taken as the constant \eqn{a} in
    the Conover-Salsburg scores.  Defaults to 5.}
  \item{j}{a numeric, the value taken as the constant \eqn{j} in the
    Koziol-Nemec scores.  Defaults to 1.}
  \item{weight}{a function whose three first arguments must correspond
    to \code{time}, \code{n.risk}, and \code{n.event} given below;
    defaults to \code{logrank_weight}.}
  \item{time}{a numeric vector, the ordered distinct time points.}
  \item{n.risk}{a numeric vector, the number of subject at risk at each
    time point specified in \code{time}.}
  \item{n.event}{a numeric vector, the number of events at each time
    point specified in \code{time}.}
  \item{type}{a character, one of \code{"logrank"} (default),
    \code{"Gehan-Breslow"}, \code{"Tarone-Ware"}, \code{"Prentice"},
    \code{"Prentice-Marek"}, \code{"Andersen-Borgan-Gill-Keiding"},
    \code{"Fleming-Harrington"} or \code{"Self"} (see
    \code{\link{surv_test}} for details).}
  \item{rho}{a numeric vector, the \eqn{\rho} constant when \code{type} is
    \code{"Tarone-Ware"}, \code{"Fleming-Harrington"} or \code{"Self"}.
    Defaults to 0.5 in the first case and 0 in the other cases.}
  \item{gamma}{a numeric vector, the \eqn{\gamma} constant when \code{type} is
    \code{"Fleming-Harrington"} or \code{"Self"}.  Defaults to 0.}
  \item{minprob}{a fraction between 0 and 0.5.}
  \item{maxprob}{a fraction between 0.5 and 1.}
  \item{scores}{a numeric vector, the scores corresponding to each level
    of an ordered factor.  Defaults to \code{1:nlevels(x)}.}
  \item{\dots}{
    \code{logrank_trafo}: further arguments to be passed to
    \code{weight}.
    % do not remove!
    \code{mcp_trafo}: factor name and contrast matrix (as matrix or
    character) in a \samp{tag = value} format for multiple comparisons
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

  \code{normal_trafo}, \code{median_trafo},\code{savage_trafo},
  \code{consal_trafo} and \code{koziol_trafo} computes normal scores,
  median scores, Savage scores, Conover-Salsburg scores (see
  \code{\link{neuropathy}}) and Koziol-Nemec scores, respectively, for
  location problems.

  \code{ansari_trafo} and \code{fligner_trafo} compute Ansari-Bradley
  and Fligner scores, respectively, for scale problems.

  \code{logrank_trafo} computes weighted logrank scores for censored
  data, allowing for user-defined weight functions through the
  \code{weight} argument (see \code{\link{GTSG}} for an example).

  A \code{trafo} function with modified default arguments is usually fed
  into \code{\link{independence_test}} via the \code{xtrafo} or
  \code{ytrafo} arguments.

  Fine tuning (different transformations for different variables) is
  possible by supplying a named list of functions to the \code{var_trafo}
  argument.

}
\note{

  In \pkg{coin} version 1.1-0 (and later), all transformation functions
  are now passing through missing values (i.e., \code{NA}s).
  Furthermore, \code{median_trafo} and \code{logrank_trafo} are now
  increasing functions of the data (in conformity with most other
  transformations in this package).

}
\value{

  A named matrix with \code{nrow(data)} rows and arbitrary number of
  columns.  User-supplied transformations must return a numeric vector
  or matrix.

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
