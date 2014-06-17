\name{ScaleTests}
\alias{taha_test}
\alias{taha_test.formula}
\alias{taha_test.IndependenceProblem}
\alias{klotz_test}
\alias{klotz_test.formula}
\alias{klotz_test.IndependenceProblem}
\alias{mood_test}
\alias{mood_test.formula}
\alias{mood_test.IndependenceProblem}
\alias{ansari_test}
\alias{ansari_test.formula}
\alias{ansari_test.IndependenceProblem}
\alias{fligner_test}
\alias{fligner_test.formula}
\alias{fligner_test.IndependenceProblem}
\alias{conover_test}
\alias{conover_test.formula}
\alias{conover_test.IndependenceProblem}

\title{ Independent Two- and K-Sample Scale Tests }

\description{
  Testing the equality of the distributions of a numeric response in two
  or more independent groups against scale alternatives.
}
\usage{

\method{taha_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{taha_test}{IndependenceProblem}(object,
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{klotz_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{klotz_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{mood_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{mood_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{ansari_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{ansari_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{fligner_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{fligner_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, \dots)

\method{conover_test}{formula}(formula, data, subset = NULL, weights = NULL, \dots)
\method{conover_test}{IndependenceProblem}(object,
    conf.int = FALSE, conf.level = 0.95, \dots)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where
    \code{y} is a numeric variable giving the data values and \code{x} a
    factor with two or more levels giving the corresponding groups.
    \code{block} is an optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}
  \item{object}{an object of class \code{"IndependenceProblem"}.}
  \item{conf.int}{a logical indicating whether a confidence interval for
    the ratio of scales should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{ties.method}{a character, two methods are available to adjust
    scores for ties, either the score generating function is applied to
    \code{"mid-ranks"} (default) or the scores computed based on random
    ranks are averaged for all tied values (\code{"average-scores"}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested.  In particular, the methods
  documented here are designed to detect scale alternatives.

  The test procedures documented here perform the Taha test
  (\code{taha_test}), the Klotz test (\code{klotz_test}), the Mood test
  (\code{mood_test}), the Ansari-Bradley test (\code{ansari_test}), the
  Fligner-Killeen test (\code{fligner_test}), and the Conover-Iman test
  (\code{conover_test}).  For a general description of these tests we
  refer to Hollander & Wolfe (1999).  The Fligner-Killeen test
  implemented here uses median centering in each of the samples, as
  suggested by Conover, Johnson & Johnson (1981), whereas the
  Conover-Iman test, following Conover and Iman (1978), uses mean
  centering in each of the samples.

  The asymptotic null distribution is computed by default for all
  procedures.  Exact \eqn{p}-values may be computed for the two-sample
  problems and can be approximated via Monte Carlo resampling for all
  procedures.  The exact \eqn{p}-values are computed either by the shift
  algorithm (Streitberg & Roehmel, 1986, 1987) or by the split-up
  algorithm (van de Wiel, 2001).

  In the two-sample case the two-sided hypothesis \eqn{var(Y_1) /
    var(Y_2) = 1} is tested, where \eqn{var(Y_i)} is the variance of the
  responses in the \eqn{i}th sample.  Confidence intervals for the ratio
  of scales are available for these tests and are computed
  according to Bauer (1972).  In case \code{alternative = "less"}, the
  null hypothesis \eqn{var(Y_1) / var(Y_2) \ge 1} is tested and
  \code{alternative = "greater"} corresponds to \eqn{var(Y_1) /
    var(Y_2) \le 1}.

  For the adjustment of scores for tied values see Hajek, Sidak and Sen
  (1999, pp. 133--135).

}

\note{
  In the two-sample case, a \emph{large} value of the Ansari-Bradley statistic
  indicates that sample 1 is \emph{less} variable than sample 2, whereas
  \emph{large} values of the statistics due to Taha, Klotz, Mood and
  Fligner-Killeen indicate that sample 1 is \emph{more} variable than
  sample 2.
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

  William J. Conover, Mark E. Johnson & Myrle M. Johnson (1981).
  A comparative study of tests for homogeneity of variances, with
  applications to the outer continental shelf bidding data.
  \emph{Technometrics} \bold{23}, 351--361.

  Conover, W. J. and Iman, R. L.  (1978).  Some Exact Tables for the
  Squared Ranks Test.  \emph{Communications in Statistics -- Simulation
    and Computation} \bold{7}(5), 491--513.

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

  ### Serum Iron Determination Using Hyland Control Sera
  ### Hollander & Wolfe (1999), page 147
  sid <- data.frame(
      serum = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
                101, 96, 97, 102, 107, 113, 116, 113, 110, 98,
                107, 108, 106, 98, 105, 103, 110, 105, 104,
                100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99),
      method = factor(gl(2, 20), labels = c("Ramsay", "Jung-Parekh")))

  ### Ansari-Bradley test, asymptotical p-value
  ansari_test(serum ~ method, data = sid)

  ### exact p-value
  ansari_test(serum ~ method, data = sid, distribution = "exact")


  ### Platelet Counts of Newborn Infants
  ### Hollander & Wolfe (1999), Table 5.4, page 171
  platalet_counts <- data.frame(
      counts = c(120, 124, 215, 90, 67, 95, 190, 180, 135, 399,
                 12, 20, 112, 32, 60, 40),
      treatment = factor(c(rep("Prednisone", 10), rep("Control", 6))))

  ### Lepage test, Hollander & Wolfe (1999), page 172
  lt <- independence_test(counts ~ treatment, data = platalet_counts,
      ytrafo = function(data) trafo(data, numeric_trafo = function(x)
          cbind(rank_trafo(x), ansari_trafo(x))),
      teststat = "quad", distribution = approximate(B = 9999))
  lt

  ### where did the rejection come from? Use maximum statistic
  ### instead of a quadratic form
  ltmax <- independence_test(counts ~ treatment, data = platalet_counts,
      ytrafo = function(data) trafo(data, numeric_trafo = function(x)
          matrix(c(rank_trafo(x), ansari_trafo(x)), ncol = 2,
                 dimnames = list(1:length(x), c("Location", "Scale")))),
      teststat = "max")

  ### single-step adjustment points to a difference in location
  pvalue(ltmax, method = "single-step")

  ### step-down is slightly more powerful
  pvalue(ltmax, method = "step-down")

  ### The same results are obtained from simple Sidak or Sidak-Holm procedures,
  ### respectively, since the correlation between the Wilcoxon and
  ### Ansari-Bradley test statistics is zero
  cov2cor(covariance(ltmax))
  pvalue(ltmax, method = "single-step", distribution = "marginal", type = "Sidak")
  pvalue(ltmax, method = "step-down", distribution = "marginal", type = "Sidak")
}
\keyword{htest}
