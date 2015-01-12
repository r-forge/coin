\name{SymmetryTests}
\alias{quade_test}
\alias{quade_test.formula}
\alias{quade_test.SymmetryProblem}
\alias{friedman_test}
\alias{friedman_test.formula}
\alias{friedman_test.SymmetryProblem}
\alias{wilcoxsign_test}
\alias{wilcoxsign_test.formula}
\alias{wilcoxsign_test.SymmetryProblem}
\alias{sign_test}
\alias{sign_test.formula}
\alias{sign_test.SymmetryProblem}

\title{ Symmetry Tests }

\description{
  Testing the symmetry of a response for repeated measurements in a
  complete block design.
}
\usage{
\method{quade_test}{formula}(formula, data, subset = NULL, \dots)
\method{quade_test}{SymmetryProblem}(object, \dots)

\method{friedman_test}{formula}(formula, data, subset = NULL, \dots)
\method{friedman_test}{SymmetryProblem}(object, \dots)

\method{wilcoxsign_test}{formula}(formula, data, subset = NULL, \dots)
\method{wilcoxsign_test}{SymmetryProblem}(object,
    zero.method = c("Pratt", "Wilcoxon"), \dots)

\method{sign_test}{formula}(formula, data, subset = NULL, \dots)
\method{sign_test}{SymmetryProblem}(object, \dots)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where
    \code{y} is a numeric variable giving the data values and \code{x}
    is a factor with two (\code{wilcoxsign_test} and \code{sign_test})
    or more levels giving the corresponding groups.  \code{block} is an
    optional factor (which is generated automatically when omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{object}{an object inheriting from class
    \code{"SymmetryProblem"}.}
  \item{zero.method}{a character specifying the way zeros, e.g., induced
    by tied repeated measurements, are handled (see \sQuote{Details}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of symmetry of \code{y} across \code{x} is
  tested.

  When \code{x} is an ordered factor \code{friedman_test} performs the
  Page test.  In this case the test statistic is a scalar and the
  alternative hypothesis can be specified using the \code{alternative}
  argument.  The default scores, \code{1:nlevels(x)}, can be altered
  using the \code{scores} argument (see \code{\link{symmetry_test}}).

  For \code{wilcoxsign_test} and \code{sign_test}, formulae of the form
  \code{y ~ x | block} and \code{y ~ x} are allowed.  The latter is
  interpreted in the sense that \code{y} is the first and \code{x} the
  second measurement on the same observation.  For
  \code{wilcoxsign_test}, the default method of handling zeros
  (\code{zero.method = "Pratt"}) due to Pratt (1959) first rank
  transforms the absolute differences (including zeros) and then
  discards the ranks corresponding to the zero differences.  The
  original proposal by Wilcoxon (\code{zero.method = "Wilcoxon"}), as
  explained in Hollander & Wolfe (1999), first discards the zero
  differences and then rank transforms the remaining absolute
  differences.

}
\note{

  Starting with \pkg{coin} version 1.0-16, the \code{zero.method}
  argument replaced the (now removed) \code{ties.method} argument.  The
  current default is \code{zero.method = "Pratt"} whereas the earlier
  versions had \code{ties.method = "HollanderWolfe"} which is now
  equivalent to \code{zero.method = "Wilcoxon"}.

}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\references{
  Quade, D.  (1979).  Using weighted rankings in the analysis of complete blocks
  with additive block effects.  \emph{Journal of the American Statistical
    Association} \bold{74}(367), 680--683.

  Myles Hollander & Douglas A. Wolfe (1999).
  \emph{Nonparametric Statistical Methods}, 2nd Edition.
  New York: John Wiley & Sons.

  John W. Pratt (1959).
  Remarks on zeros and ties in the Wilcoxon signed rank procedures.
  \emph{Journal of the American Statistical Association} \bold{54}(287),
  655--667.

  Bernd Streitberg & Joachim Roehmel (1986).
  Exact distributions for permutations and rank tests: An introduction
  to some recently published algorithms.
  \emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

  Bernd Streitberg & Joachim Roehmel (1987).
  Exakte Verteilungen fuer Rang- und Randomisierungstests im allgemeinen
  c-Stichprobenfall.
  \emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

}
\examples{
  ### Data from Quade (1979, p. 683)
  dta <- data.frame(
      y = c(52, 45, 38,
            63, 79, 50,
            45, 57, 39,
            53, 51, 43,
            47, 50, 56,
            62, 72, 49,
            49, 52, 40),
       x = factor(rep(LETTERS[1:3], 7)),
       b = factor(rep(1:7, each = 3)))

  ### Approximate (Monte Carlo) Friedman test
  friedman_test(y ~ x | b, data = dta,
                distribution = approximate(B = 10000))

  ### Approximate (Monte Carlo) Quade test
  ### Quade (1979, p.683), W = 8.157
  (qt <- quade_test(y ~ x | b, data = dta,
                    distribution = approximate(B = 10000)))

  ### Comparison with R's quade.test() function
  quade.test(y ~ x | b, data = dta)

  ### quade.test() uses an F-statistic
  b <- nlevels(qt@statistic@block)
  A <- sum(qt@statistic@y^2)
  B <- sum(statistic(qt, "linear")^2) / b
  (b - 1) * B / (A - B) # F = 8.3765


  ### Hollander & Wolfe (1999), Table 7.1, page 274
  ### Comparison of three methods ("round out", "narrow angle", and
  ###  "wide angle") for rounding first base.
  RoundingTimes <- data.frame(
      times = c(5.40, 5.50, 5.55,
                5.85, 5.70, 5.75,
                5.20, 5.60, 5.50,
                5.55, 5.50, 5.40,
                5.90, 5.85, 5.70,
                5.45, 5.55, 5.60,
                5.40, 5.40, 5.35,
                5.45, 5.50, 5.35,
                5.25, 5.15, 5.00,
                5.85, 5.80, 5.70,
                5.25, 5.20, 5.10,
                5.65, 5.55, 5.45,
                5.60, 5.35, 5.45,
                5.05, 5.00, 4.95,
                5.50, 5.50, 5.40,
                5.45, 5.55, 5.50,
                5.55, 5.55, 5.35,
                5.45, 5.50, 5.55,
                5.50, 5.45, 5.25,
                5.65, 5.60, 5.40,
                5.70, 5.65, 5.55,
                6.30, 6.30, 6.25),
      methods = factor(rep(c("Round Out", "Narrow Angle", "Wide Angle"), 22)),
      block = factor(rep(1:22, rep(3, 22))))

  ### classical global test
  friedman_test(times ~ methods | block, data = RoundingTimes)

  ### parallel coordinates plot
  matplot(t(matrix(RoundingTimes$times, ncol = 3, byrow = TRUE)),
          type = "l", col = 1, lty = 1, axes = FALSE, ylab = "Time",
          xlim = c(0.5, 3.5))
  axis(1, at = 1:3, labels = c("Round Out", "Narrow Angle", "Wide Angle"))
  axis(2)

  ### where do the differences come from?
  ### Wilcoxon-Nemenyi-McDonald-Thompson test
  ### Hollander & Wolfe (1999), page 295
  ### all pairwise comparisons
  rtt <- symmetry_test(times ~ methods | block, data = RoundingTimes,
       xtrafo = mcp_trafo(methods = "Tukey"),
       ytrafo = function(data)
           trafo(data, numeric_trafo = rank_trafo, block = RoundingTimes$block)
       )

  ### a global test, again
  print(pvalue(rtt))

  ### simultaneous P-values for all pair comparisons
  ### Wide Angle vs. Round Out differ (Hollander and Wolfe, 1999, page 296)
  print(pvalue(rtt, method = "single-step"))


  ### Strength Index of Cotton, Hollander & Wolfe (1999), Table 7.5, page 286
  sc <- data.frame(block = factor(c(rep(1, 5), rep(2, 5), rep(3, 5))),
                   potash = ordered(rep(c(144, 108, 72, 54, 36), 3),
                                    levels = c(144, 108, 72, 54, 36)),
                   strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                                7.68, 7.57, 7.73, 8.15, 8.00,
                                7.21, 7.80, 7.74, 7.87, 7.93))

  ### Page test for ordered alternatives
  friedman_test(strength ~ potash | block, data = sc)

  ### one-sided p-value
  pvalue(friedman_test(strength ~ potash | block, data = sc,
                       alternative = "greater"))

  ### approximate null distribution via Monte Carlo
  pvalue(friedman_test(strength ~ potash | block, data = sc,
                       distribution = approximate(B = 9999)))


  ### example from ?wilcox.test
  y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
  y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

  wilcox.test(y1, y2, paired = TRUE, alternative = "greater")
  (wt <- wilcoxsign_test(y1 ~ y2, alternative = "greater",
                         distribution = "exact"))
  statistic(wt, "linear") # same as wilcox.test
  midpvalue(wt) # mid-p-value

  (st <- sign_test(y1 ~ y2, alternative = "greater",
                   distribution = "exact"))
  midpvalue(st) # mid-p-value

  ### with explicit group and block information
  dta <- data.frame(y = c(y1, y2), x = gl(2, length(y1)),
                    block = factor(rep(seq_along(y1), 2)))

  wilcoxsign_test(y ~ x | block, data = dta,
                  alternative = "greater", distribution = "exact")

  sign_test(y ~ x | block, data = dta,
            alternative = "greater", distribution = "exact")

  ### for two samples, the Quade test is equivalent to the signed-rank test
  quade_test(y ~ x | block, data = dta, distribution = "exact")

  wilcoxsign_test(y ~ x | block, data = dta, distribution = "exact")

  ### for two samples, the Friedman test is equivalent to the sign test
  friedman_test(y ~ x | block, data = dta, distribution = "exact")

  sign_test(y ~ x | block, data = dta, distribution = "exact")
}
\keyword{htest}
