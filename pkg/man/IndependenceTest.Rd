\name{IndependenceTest}
\alias{independence_test}
\alias{independence_test.formula}
\alias{independence_test.table}
\alias{independence_test.IndependenceProblem}

\title{ General Independence Tests }

\description{
  The independence between two sets of variables of arbitrary
  measurement scales, possibly stratified in blocks, is tested
  conditional on the data.
}
\usage{
\method{independence_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{independence_test}{table}(object, \dots)
\method{independence_test}{IndependenceProblem}(object,
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL,
    check = NULL, paired = FALSE, \dots)
}
\arguments{
  \item{formula}{a formula of the form
    \code{y1 + ... + yp ~ x1 + ... + xq | block} where the variables on
    the left and right hand side may be measured on arbitrary scales
    (including censored ones on the left hand side) and \code{block} is
    an optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.  Alternatively, an object of class
    \code{\link[Biobase]{class.ExpressionSet}} may be specified.  In
    this case, all variables in \code{formula}, except \code{.}, are
    first evaluated in the \code{pData} data frame.  The dot (\code{.})
    refers to the matrix of expression levels (\code{exprs} slot).}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}
  \item{object}{an object inheriting from class
    \code{"IndependenceProblem"} or an object of class \code{"table"}.}
  \item{teststat}{a character, the type of test statistic to be applied:
    either a standardized scalar test statistic (\code{"scalar"}), or a
    maximum type statistic (\code{"max"}, default) or a quadratic form
    (\code{"quad"}).}
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
    \code{"NullDistribution"}.}
  \item{alternative}{a character, the alternative hypothesis must be one
    of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.}
  \item{xtrafo}{a function of transformations (see \code{\link{trafo}})
    to be applied to the variables on the right hand side of
    \code{formula}, see below.}
  \item{ytrafo}{a function of transformations (see \code{\link{trafo}})
    to be applied to the variables on the left hand side of
    \code{formula}, see below.}
  \item{scores}{a named list of scores to be attached to ordered
    factors. In case a variable is an unordered factor, it is coerced to
    \code{ordered} first.}
  \item{check}{a function to be applied to objects of class
    \code{"IndependenceTest"} in order to check for specific properties
    of the data.}
  \item{paired}{a logical, indicating that paired data have been
    transformed in such a way that the resulting unstandardized linear
    statistic is the sum of the absolute values of the positive
    differences between the paired observations.}
  \item{\dots}{further arguments to be passed to or from
    methods. Currently, none of the additional arguments is passed to
    any function.}
}
\details{

  The null hypothesis of the independence between the variables on the
  left hand side and the variables on the right hand side of
  \code{formula}, possibly stratified by \code{block}, is tested.  The
  vector supplied via the \code{weights} argument is interpreted as
  observation counts.

  This function is the basic workhorse called by all other convenience
  functions, mainly by supplying transformations via the \code{xtrafo}
  argument and influence functions via the \code{ytrafo} argument.

  The \code{scores} argument leads to linear-by-linear association tests
  against ordered alternatives.  If the formula \code{y ~ x} was
  supplied and both \code{y} and \code{x} are factors,
  \code{scores = list(y = 1:k, x = c(1, 4, 6))} first triggers a
  coercion to class \code{ordered} of both variables and attaches the
  list elements as scores.  The length of a score vector needs to be
  equal the number of levels of the factor of interest.

  The basis of this function is the framework for conditional inference
  procedures by Strasser & Weber (1999).  The theory and this
  implementation are explained and illustrated in Hothorn, Hornik, van
  de Wiel & Zeileis (2006).

}
\note{
  Starting with \pkg{coin} version 1.1-0, maximum-type statistics and
  quadratic forms can no longer be specified using \code{teststat =
    "maxtype"} and \code{teststat = "quadtype"} respectively (as was
  used in versions prior to 0.4-5).
}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\references{

  Helmut Strasser & Christian Weber (1999).
  On the asymptotic theory of permutation statistics.
  \emph{Mathematical Methods of Statistics} \bold{8}(2), 220--250.

  Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel & Achim Zeileis
  (2006).
  A Lego System for Conditional Inference.
  \emph{The American Statistician} \bold{60}(3), 257--263.

  Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel & Achim Zeileis
  (2008).
  Implementing a class of permutation tests: The coin package.
  \emph{Journal of Statistical Software} \bold{28}(8), 1--23.
  \url{http://www.jstatsoft.org/v28/i08/}

}
\examples{

  ### independence of asat and group via normal scores test
  independence_test(asat ~ group, data = asat,

    ### exact null distribution
    distribution = "exact",

    ### one-sided test
    alternative = "greater",

    ### apply normal scores to asat$asat
    ytrafo = function(data) trafo(data, numeric_trafo = normal_trafo),

    ### indicator matrix of 1st level of group
    xtrafo = function(data) trafo(data, factor_trafo = function(x)
        matrix(x == levels(x)[1], ncol = 1))
  )

  ### same as
  normal_test(asat ~ group, data = asat, distribution = "exact",
              alternative = "greater")


  ### permutation test for paired observations
  y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
  y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
  dta <- data.frame(y = c(y1, y2), x = gl(2, length(y1)),
                    block = factor(rep(seq_along(y1), 2)))
  independence_test(y ~ x | block, data = dta,
                    distribution = "exact", alternative = "greater")

  ### alternatively: transform data and set 'paired = TRUE'
  diff <- y1 - y2
  y <- as.vector(rbind(abs(diff) * (diff >= 0), abs(diff) * (diff < 0)))
  x <- factor(rep(0:1, length(diff)), labels = c("pos", "neg"))
  block <- gl(length(diff), 2)
  independence_test(y ~ x | block,
                    distribution = "exact", alternative = "greater",
                    paired = TRUE)


  ### if you are interested in the internals:
  \dontrun{
      browseURL(system.file("documentation", "html", "index.html",
                            package = "coin"))
  }

}
\keyword{htest}
