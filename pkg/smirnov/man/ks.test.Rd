% File src/library/stats/man/ks.test.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2020 R Core Team
% Distributed under GPL 2 or later

\name{ks.test}
\alias{ks.test}
\alias{ks.test.default}
\alias{ks.test.formula}
\alias{psmirnov}
\alias{qsmirnov}
\title{Kolmogorov-Smirnov and Smirnov Tests}
\description{
  Perform a Kolmogorov-Smirnov (one-sample) or Smirnov (two-sample) test.
}
\usage{
ks.test(x, \dots)
\method{ks.test}{default}(x, y, \dots,
        alternative = c("two.sided", "less", "greater"),
        exact = NULL, simulate.p.value = FALSE, B = 2000)
\method{ks.test}{formula}(formula, data, subset, na.action, \dots)
psmirnov(q, m, n = length(z) - m, z = NULL, 
         two.sided = TRUE, exact = TRUE, simulate = FALSE, B = 2000,
         lower.tail = TRUE, log.p = FALSE)
qsmirnov(p, m, n = length(z) - m, z = NULL, two.sided = TRUE, \dots)
}
\arguments{
  \item{x}{a numeric vector of data values.}
  \item{y}{either a numeric vector of data values, or a character string
    naming a cumulative distribution function or an actual cumulative
    distribution function such as \code{pnorm}.  Only continuous CDFs
    are valid.}
  \item{\dots}{parameters of the distribution specified (as a character
    string) by \code{y} for \code{ks.test} or arguments of \code{psmirnov}
    for \code{qsmirnov}.}
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"two.sided"} (default), \code{"less"}, or
    \code{"greater"}.  You can specify just the initial letter of the
    value, but the argument name must be given in full.
    See \sQuote{Details} for the meanings of the possible values.}
  \item{exact}{\code{NULL} or a logical indicating whether an exact
    p-value should be computed.  See \sQuote{Details} for the meaning of
    \code{NULL}.}
  \item{simulate.p.value,simulate}{a logical indicating whether to compute p-values by
          Monte Carlo simulation.}
  \item{B}{an integer specifying the number of replicates used in the
          Monte Carlo test.}
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs}
    is a numeric variable giving the data values and \code{rhs} either
    \code{1} for a one-sample Kolmogorov-Smirnov test or a factor
    with two levels giving the corresponding groups for a two-sample Smirnov
    test.}
  \item{data}{an optional matrix or data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{q}{a numeric vector of quantiles.}
  \item{p}{a numeric vector of probabilities.}
  \item{m}{length of \code{x}-vector, sample size in the first group in
    the two-sample case.}
  \item{n}{length of \code{y}-vector, sample size in the second group in
    the two-sample case.}
  \item{z}{a numeric vector of all data values (\code{x} and \code{y})
             when the exact conditional distribution of the Smirnov test
             statistic given the data shall be computed.}
  \item{two.sided}{a logical indicating whether absolute (\code{TRUE}) or
                   raw differences of frequencies define the test statistic.}
  \item{lower.tail}{a logical, if \code{TRUE} (default), probabilities are
                    \eqn{P[D < q]}, otherwise, \eqn{P[D \ge q]}.}
  \item{log.p}{a logical, if \code{TRUE} (default), probabilities are given
               as log-probabilities.}

}
\details{
  If \code{y} is numeric, a two-sample Smirnov test of the null hypothesis
  that \code{x} and \code{y} were drawn from the same distribution is performed. 

  Alternatively, \code{y} can be a character string naming a continuous
  (cumulative) distribution function, or such a function.  In this case,
  a one-sample Kolmogorov-Smirnov test is carried out of the null that the distribution
  function which generated \code{x} is distribution \code{y} with
  parameters specified by \code{\dots}. The presence of ties always generates a 
  warning in the one-sided case, since continuous
  distributions do not generate them.  If the ties arose from rounding
  the tests may be approximately valid, but even modest amounts of
  rounding can have a significant effect on the calculated statistic.

  Missing values are silently omitted from \code{x} and (in the
  two-sample case) \code{y}.

  The possible values \code{"two.sided"}, \code{"less"} and
  \code{"greater"} of \code{alternative} specify the null hypothesis
  that the true distribution function of \code{x} is equal to, not less
  than or not greater than the hypothesized distribution function
  (one-sample case) or the distribution function of \code{y} (two-sample
  case), respectively.  This is a comparison of cumulative distribution
  functions, and the test statistic is the maximum difference in value,
  with the statistic in the \code{"greater"} alternative being
  \eqn{D^+ = \max_u [ F_x(u) - F_y(u) ]}{D^+ = max[F_x(u) - F_y(u)]}.
  Thus in the two-sample case \code{alternative = "greater"} includes
  distributions for which \code{x} is stochastically \emph{smaller} than
  \code{y} (the CDF of \code{x} lies above and hence to the left of that
  for \code{y}), in contrast to \code{\link{t.test}} or
  \code{\link{wilcox.test}}.

  Exact p-values are not available for the one-sided case in the presence of ties.  
  If \code{exact = NULL} (the default), an
  exact p-value is computed if the sample size is less than 100 in the
  one-sample case \emph{and there are no ties}, and if the product of
  the sample sizes is less than 10000 in the two-sample case, with or
  without ties (using the algorithm described in Schröer and Trenkler, 1995).
  Otherwise, asymptotic distributions are used whose approximations may
  be inaccurate in small samples.  In the one-sample two-sided case,
  exact p-values are obtained as described in Marsaglia, Tsang & Wang
  (2003) (but not using the optional approximation in the right tail, so
  this can be slow for small p-values).  The formula of Birnbaum &
  Tingey (1951) is used for the one-sample one-sided case.

  If a one-sample Kolmogorov-Smirnov test is used, the parameters specified in
  \code{\dots} must be pre-specified and not estimated from the data.
  There is some more refined distribution theory for the KS test with
  estimated parameters (see Durbin, 1973), but that is not implemented
  in \code{ks.test}.

  \code{psmirnov} and \code{qsmirnov} compute the distribution and quantile
  function of the Smirnov test for two samples.

}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\source{
  The two-sided one-sample distribution comes \emph{via}
  Marsaglia, Tsang and Wang (2003).

  Exact distributions for the two-sample Smirnov test are computed
  by the algorithm proposed by Schröer (1991) and Schröer and Trenkler
  (1995).

}
\references{

  Vance W. Berger and YanYan Zhou  (2014).  
  Kolmogorov–Smirnov Test: Overview.  In Wiley StatsRef: \emph{Statistics Reference
  Online} (eds N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch, F. Ruggeri 
  and J.L. Teugels). \doi{10.1002/9781118445112.stat06558}.

  Zygmunt W. Birnbaum and Fred H. Tingey (1951).
  One-sided confidence contours for probability distribution functions.
  \emph{The Annals of Mathematical Statistics}, \bold{22}/4, 592--596.
  \doi{10.1214/aoms/1177729550}.

  William J. Conover (1971).
  \emph{Practical Nonparametric Statistics}.
  New York: John Wiley & Sons.
  Pages 295--301 (one-sample Kolmogorov test),
  309--314 (two-sample Smirnov test).

  James Durbin (1973).
  \emph{Distribution theory for tests based on the sample distribution
    function}.
  SIAM.

  George Marsaglia, Wai Wan Tsang and Jingbo Wang (2003).
  Evaluating Kolmogorov's distribution.
  \emph{Journal of Statistical Software}, \bold{8}/18.
  \doi{10.18637/jss.v008.i18}.

  Gunar Schröer and Dietrich Trenkler (1995). Exact and Randomization 
  Distributions of Kolmogorov-Smirnov Tests for Two or Three
  Samples. \emph{Computational Statistics & Data Analysis}, \bold{20}(2),
  185--202. \doi{10.1016/0167-9473(94)00040-P}.

}
\seealso{
  \code{\link{shapiro.test}} which performs the Shapiro-Wilk test for
  normality.
}
\examples{
require("graphics")

x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.test(x, y)
# Does x come from a shifted gamma distribution with shape 3 and rate 2?
ks.test(x+2, "pgamma", 3, 2) # two-sided, exact
ks.test(x+2, "pgamma", 3, 2, exact = FALSE)
ks.test(x+2, "pgamma", 3, 2, alternative = "gr")

# test if x is stochastically larger than x2
x2 <- rnorm(50, -1)
plot(ecdf(x), xlim = range(c(x, x2)))
plot(ecdf(x2), add = TRUE, lty = "dashed")
t.test(x, x2, alternative = "g")
wilcox.test(x, x2, alternative = "g")
ks.test(x, x2, alternative = "l")

# with ties, example from Schröer and Trenkler (1995)
# D = 3 / 7, p = 0.2424242
ks.test(c(1, 2, 2, 3, 3), c(1, 2, 3, 3, 4, 5, 6), exact = TRUE)

# formula interface, see ?wilcox.test
kst <- ks.test(Ozone ~ Month, data = airquality,
               subset = Month \%in\% c(5, 8))
# quantile-quantile plot + confidence bands
plot(confband(kst)) # => null hypothesis not plausible, 
                    # shift alternative not plausible
}
\keyword{htest}
