\name{LocationTests}
\alias{oneway_test}
\alias{oneway_test.formula}
\alias{oneway_test.IndependenceProblem}
\alias{wilcox_test}
\alias{wilcox_test.formula}
\alias{wilcox_test.IndependenceProblem}
\alias{kruskal_test}
\alias{kruskal_test.formula}
\alias{kruskal_test.IndependenceProblem}
\alias{normal_test}
\alias{normal_test.formula}
\alias{normal_test.IndependenceProblem}
\alias{median_test}
\alias{median_test.formula}
\alias{median_test.IndependenceProblem}
\alias{savage_test}
\alias{savage_test.formula}
\alias{savage_test.IndependenceProblem}

\title{ Independent Two- and K-Sample Location Tests }

\description{
  Testing the equality of the distributions of a numeric response in two
  or more independent groups against shift alternatives.
}
\usage{

\method{oneway_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{oneway_test}{IndependenceProblem}(object, \dots)

\method{wilcox_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{wilcox_test}{IndependenceProblem}(object,
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{kruskal_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{kruskal_test}{IndependenceProblem}(object,
    distribution = c("asymptotic", "approximate"), \dots)

\method{normal_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{normal_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{median_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{median_test}{IndependenceProblem}(object,
    mid.score = c("0", "0.5", "1"),
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{savage_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{savage_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where
    \code{y} is a numeric variable giving the data values and \code{x} a
    factor with two or more levels giving the corresponding
    groups. \code{block} is an optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}
  \item{object}{an object of class \code{IndependenceProblem}.}
  \item{conf.int}{a logical indicating whether a confidence interval for
    the difference in location should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{distribution}{a character, the null distribution of the test
    statistic can be computed exactly (\code{"exact"}) or can be
    approximated by its asymptotic distribution (\code{"asymptotic"}) or
    via Monte-Carlo resampling (\code{"approximate"}). Alternatively,
    the functions \code{\link{exact}}, \code{\link{approximate}} or
    \code{\link{asymptotic}} can be used to specify how the exact
    conditional distribution of the test statistic should be calculated
    or approximated.}
  \item{ties.method}{a character, two methods are available to adjust
    scores for ties: either the score generating function is applied to
    \code{"mid-ranks"} or the scores computed based on random ranks are
    averaged for all tied values (\code{"average-scores"}).}
  \item{mid.score}{a character, the score assigned to observations
    exactly equal to the median.  Must be one of \code{"0"} (default),
    \code{"0.5"} or \code{"1"} (see \sQuote{Details}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested. In particular, the methods
  documented here are designed to detect shift alternatives. For a
  general description of the test procedures documented here we refer to
  Hollander & Wolfe (1999).

  The test procedures apply a rank-based transformation, i.e., either
  the rank transformation itself or a function of the ranks, to the
  response values \code{y}, except \code{oneway_test} which computes its
  test statistic using the untransformed response values.

  The asymptotic null distribution is computed by default for all
  procedures.  Exact \eqn{p}-values may be computed for the two-sample
  problems and can be approximated via Monte Carlo resampling for all
  procedures.  The exact \eqn{p}-values are computed either by the shift
  algorithm (Streitberg & Roehmel, 1986, 1987) or by the split-up
  algorithm (van de Wiel, 2001).

  In the two-sample case the linear rank tests \code{wilcox_test},
  \code{normal_test}, \code{median_test} and \code{savage_test} can be
  used to test the two-sided hypothesis \eqn{H_0: \mu = 0}, where
  \eqn{\mu = Y_1 - Y_2} and \eqn{Y_i} is the median of the responses in
  the \eqn{i}th sample.  Confidence intervals for the difference in
  location are available for these rank-based procedures and are
  computed according to Bauer (1972).  In case
  \code{alternative = "less"}, the null hypothesis \eqn{H_0: \mu \ge 0}
  is tested and \code{alternative = "greater"} corresponds to the null
  hypothesis \eqn{H_0: \mu \le 0}.

  The median test offers a choice of mid-score, i.e., the score assigned
  to observations exactly equal to the median, using \code{"mid-score"}.
  In the two-sample case \code{"0"} implies that the linear test
  statistic is simply the number of subjects in the 2nd sample with
  observations greater than the median of the pooled sample.  The linear
  test statistic for the last alternative \code{"1"} is the number of
  subjects in the 2nd sample whose observations are greater or equal to
  the median of the pooled sample.  If \code{"0.5"} is selected, the
  linear test statistic is the mean of the test statistics corresponding
  to the first and last alternatives and has a symmetric distribution,
  or at least approximately so, under the null hypothesis (see Hajek,
  Sidak and Sen, 1999, pp. 97--98).

  In case \code{x} is an ordered factor, the linear-by-linear
  association test for ordered alternatives is computed.

  For the adjustment of scores for tied values see Hajek, Sidak and Sen
  (1999, pp. 133--135).

}
\value{

  An object inheriting from class \code{\link{IndependenceTest-class}} with
  methods \code{\link{show}}, \code{\link{statistic}}, \code{\link{expectation}},
  \code{\link{covariance}} and \code{\link{pvalue}}. The null distribution
  can be inspected by \code{\link{pperm}}, \code{\link{dperm}},
  \code{\link{qperm}} and \code{\link{support}} methods. Confidence
  intervals can be extracted by \code{confint}.

}
\references{

  Myles Hollander & Douglas A. Wolfe (1999).
  \emph{Nonparametric Statistical Methods}, 2nd Edition.
  New York: John Wiley & Sons.

  Bernd Streitberg & Joachim Roehmel (1986).
  Exact distributions for permutations and rank tests: An introduction
  to some recently published algorithms.
  \emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

  Bernd Streitberg & Joachim Roehmel (1987).
  Exakte Verteilungen fuer Rang- und Randomisierungstests im allgemeinen
  c-Stichprobenfall.
  \emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

  Mark A. van de Wiel (2001).
  The split-up algorithm: a fast symbolic method for computing p-values
  of rank statistics.
  \emph{Computational Statistics} \bold{16}, 519--538.

  David F. Bauer (1972).
  Constructing confidence sets using rank statistics.
  \emph{Journal of the American Statistical Association} \bold{67},
  687--690.

  Jaroslav Hajek, Zbynek Sidak & Pranab K. Sen (1999),
  \emph{Theory of Rank Tests}, 2nd Edition.
  San Diego: Academic Press.

}
\examples{

  ### Tritiated Water Diffusion Across Human Chorioamnion
  ### Hollander & Wolfe (1999), Table 4.1, page 110
  water_transfer <- data.frame(
      pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
             1.15, 0.88, 0.90, 0.74, 1.21),
      age = factor(c(rep("At term", 10), rep("12-26 Weeks", 5))))

  ### Wilcoxon-Mann-Whitney test, cf. Hollander & Wolfe (1999), page 111
  ### exact p-value and confidence interval for the difference in location
  ### (At term - 12-26 Weeks)
  wt <- wilcox_test(pd ~ age, data = water_transfer,
                    distribution = "exact", conf.int = TRUE)
  print(wt)

  ### extract observed Wilcoxon statistic, i.e, the sum of the
  ### ranks for age = "12-26 Weeks"
  statistic(wt, "linear")

  ### its expectation
  expectation(wt)

  ### and variance
  covariance(wt)

  ### and the exact two-sided p-value
  pvalue(wt)

  ##d and, finally, the confidence interval
  confint(wt)

  ### Confidence interval for difference (12-26 Weeks - At term)
  wilcox_test(pd ~ age, data = water_transfer,
              xtrafo = function(data)
                  trafo(data, factor_trafo = function(x)
                      as.numeric(x == levels(x)[2])),
              distribution = "exact", conf.int = TRUE)

  ### Permutation test, asymptotic p-value
  oneway_test(pd ~ age, data = water_transfer)

  ### approximate p-value (with 99\% confidence interval)
  pvalue(oneway_test(pd ~ age, data = water_transfer,
                     distribution = approximate(B = 9999)))
  ### exact p-value
  pt <- oneway_test(pd ~ age, data = water_transfer, distribution = "exact")
  pvalue(pt)

  ### plot density and distribution of the standardized test statistic
  layout(matrix(1:2, nrow = 2))
  s <- support(pt)
  d <- sapply(s, function(x) dperm(pt, x))
  p <- sapply(s, function(x) pperm(pt, x))
  plot(s, d, type = "S", xlab = "Test Statistic", ylab = "Density")
  plot(s, p, type = "S", xlab = "Test Statistic", ylab = "Cum. Probability")
  layout(1)


  ### Median test with different mid-scores
  ex <- data.frame(y = c(3, 4, 8, 9, 1, 2, 5, 6, 7),
                   x = factor(rep(c("no", "yes"), c(4, 5))))

  boxplot(y ~ x, data = ex)

  (mt1 <- median_test(y ~ x, data = ex, distribution = "exact"))
  (mt2 <- median_test(y ~ x, data = ex, distribution = "exact",
                      mid.score = "0.5"))
  (mt3 <- median_test(y ~ x, data = ex, distribution = "exact",
                      mid.score = "1")) # sign change!

  ### plot density and distribution of the standardized test statistics
  layout(matrix(1:3, nrow = 3))
  s1 <- support(mt1); d1 <- dperm(mt1, s1)
  plot(s1, d1, type = "h", main = "Mid-score: 0",
       xlab = "Test Statistic", ylab = "Density")
  s2 <- support(mt2); d2 <- dperm(mt2, s2)
  plot(s2, d2, type = "h", main = "Mid-score: 0.5",
       xlab = "Test Statistic", ylab = "Density")
  s3 <- support(mt3); d3 <- dperm(mt3, s3)
  plot(s3, d3, type = "h", main = "Mid-score: 1",
       xlab = "Test Statistic", ylab = "Density")
  layout(1)


  ### Length of YOY Gizzard Shad from Kokosing Lake, Ohio,
  ### sampled in Summer 1984, Hollander & Wolfe (1999), Table 6.3, page 200
  YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
                               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
                               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
                               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                    site = factor(c(rep("I", 10), rep("II", 10),
                                    rep("III", 10), rep("IV", 10))))

  ### Kruskal-Wallis test, approximate exact p-value
  kw <- kruskal_test(length ~ site, data = YOY,
                     distribution = approximate(B = 9999))
  kw
  pvalue(kw)

  ### Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
  ### Hollander & Wolfe (1999), page 244
  ### (where Steel-Dwass results are given)
  NDWD <- oneway_test(length ~ site, data = YOY,
      ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
      xtrafo = function(data) trafo(data, factor_trafo = function(x)
          model.matrix(~x - 1) \%*\% t(contrMat(table(x), "Tukey"))),
      teststat = "max", distribution = approximate(B = 90000))

  ### global p-value
  print(pvalue(NDWD))

  ### sites (I = II) != (III = IV) at alpha = 0.01 (page 244)
  print(pvalue(NDWD, method = "single-step"))


  ### Jonckheere-Terpstra test for ordered groups
  control <- c(40, 35, 38, 43, 44, 41)
  rough <- c(38, 40, 47, 44, 40, 42)
  accurate <- c(48, 40, 45, 43, 46, 44)

  pieces <- list(control, rough, accurate)
  n <- c(6, 6, 6)
  grp <- as.ordered(factor(rep(1:length(n), n)))

  (y <- unlist(pieces))
  k <- length(pieces)
  (x <- as.ordered(factor(rep(1:k, n))))

  ### look at K. The second line just sums up.
  ff <- function(x) {
      K <- contrMat(table(x), "Tukey")[, x]
      as.vector(rep(1, nrow(K)) \%*\% K)
  }

  it <- independence_test(y ~ x,
      ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo),
      xtrafo = function(data) trafo(data, ordered_trafo = ff),
      alternative = "greater")
  it

}
\keyword{htest}
