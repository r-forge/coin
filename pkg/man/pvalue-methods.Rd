\name{pvalue-methods}
\docType{methods}
\alias{pvalue}
\alias{pvalue-methods}
\alias{pvalue,NullDistribution-method}
\alias{pvalue,IndependenceTest-method}
\alias{pvalue,ScalarIndependenceTest-method}
\alias{pvalue,MaxTypeIndependenceTest-method}
\alias{pvalue,QuadTypeIndependenceTest-method}
\title{ Extract \eqn{p}-Values }

\encoding{UTF-8}

\description{

  Extracts the \eqn{p}-value from objects representing null
  distributions of independence tests.

}
\usage{
pvalue(object, \dots)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{\dots}{additional arguments:

    \code{method}: a character string, the method used for the
    \eqn{p}-value computation. Must be one of \code{"global"} (default),
    \code{"single-step"}, \code{"step-down"}, \code{"npmcp"} or
    \code{"unadjusted"}.

    \code{distribution}: a character string, the distribution
    used for the computation of adjusted \eqn{p}-values. Must be one of
    \code{"joint"} (default) or \code{"marginal"}.

    \code{type}: a character string, the type of adjustment when the
    marginal distributions are used. Must be one of \code{"Bonferroni"}
    (default) or \code{"Sidak"}.
  }
}
\section{Methods}{
  \describe{
    \item{pvalue}{extracts the \eqn{p}-value from the specified
      object.}
  }
}
\details{

  The global \eqn{p}-value ({\code{method = "global"}}) is returned by
  default, and is given with an associated 99 \% confidence interval when
  resampling is used to determine the null distribution (which for
  maximum-type statistics may be true even in the asymptotic case).

  Assuming \emph{subset pivotality}, single-step or \emph{free}
  step-down adjusted \eqn{p}-values using max-\eqn{T} procedures are
  obtained by setting \code{method} to \code{"single-step"} or
  \code{"step-down"} respectively. In both cases, the
  \code{distribution} argument specifies whether the adjustment is based
  on the joint distribution (\code{"joint"}) or the marginal
  distributions (\code{"marginal"}) of the test statistics. For
  procedures based on the marginal distributions, Bonferroni- or
  \enc{Šidák}{Sidak}-type adjustment can be specified through the
  \code{type} argument by setting it to \code{"Bonferroni"} or
  \code{"Sidak"} respectively.

  For multiple comparisons of levels of a single factor (defined via
  \code{\link{mcp_trafo}}), adjusted p-values, less affected by failure
  of the \emph{subset pivotality} condition, can be computed using a
  joint distribution-based \emph{restricted} step-down max-\eqn{T}
  procedure via \code{method = "npmcp"}. This feature is \strong{highly}
  experimental.

  The \eqn{p}-value adjustment procedures based on the joint
  distribution of the test statistics fully utilizes distributional
  characteristics, such as discreteness and dependence structure,
  whereas procedures using the marginal distributions only incorporate
  discreteness. Hence, the joint distribution-based procedures are
  typically more powerful. Details regarding the single-step and
  \emph{free} step-down procedures based on the joint distribution can
  be found in Westfall & Young (1993): in particular, this
  implementation uses Equation 2.8 with Algorithm 2.5 and 2.8
  respectively. The joint distribution \emph{restricted} step-down
  procedure is based on Westfall (1997). Westfall & Wolfinger (1997)
  provide details of the marginal distributions-based single-step and
  \emph{free} step-down procedures. The generalization of Westfall &
  Wolfinger (1997) to arbitrary test statistics, as implemented here, is
  given by Westfall & Troendle (2008)

  Unadjusted \eqn{p}-values are obtained using \code{method = "unadjusted"}.

}

\note{

  The familywise error rate (FWER) is always controlled under the global
  null hypothesis, i.e. in the \emph{weak} sense, implying that the
  smallest adjusted \eqn{p}-value is valid without further
  assumptions. Control of the FWER under any partial configuration of
  the null hypotheses, i.e. in the \emph{strong} sense, as is typically
  desired for multiple tests and comparisons, requires that the
  \emph{subset pivotality} condition holds (Westfall & Young, 1993,
  pp. 42--43; Bretz, Hothorn & Westfall, 2011, pp. 136--137). In
  addition, for methods based on the joint distribution of the test
  statistics, failure of the \emph{joint exchangeability} assumption
  (Westfall & Troendle, 2008; Bretz, Hothorn & Westfall, 2011,
  pp. 129--130) may cause excess Type I errors.

  In versions of \pkg{coin} prior to 1.1-0, a min-\eqn{P} procedure
  computing \enc{Šidák}{Sidak} single-step adjusted \eqn{p}-values
  accounting for discreteness was available when specifying
  \code{method = "discrete"}. \strong{This is now deprecated and will be
  removed in a future release} due to the introduction of a more general
  max-\eqn{T} version of the same algorithm.

}
\references{

  Frank Bretz, Torsten Hothorn & Peter Westfall (2011).
  \emph{Multiple Comparisons Using R}.
  Boca Raton: CRC Press.

  Peter H. Westfall (1997).
  Multiple Testing of General Contrasts Using Logical Constraints and
  Correlations.
  \emph{Journal of the American Statistical Association} \bold{92}, 299--306.

  Peter H. Westfall & James F. Troendle (2008).
  Multiple Testing with Minimal Assumptions.
  \emph{Biometrical Journal} \bold{50}, 745--755.

  Peter H. Westfall & Russell D. Wolfinger (1997).
  Multiple Tests with Discrete Distributions.
  \emph{The American Statistician} \bold{51}, 3--8.

  Peter H. Westfall & Stanley S. Young (1993).
  \emph{Resampling-based Multiple Testing: Examples and Methods for
  \eqn{p}-Value Adjustment}.
  New York: John Wiley and Sons.

}
\examples{

  ### artificial 2-sample problem
  df <- data.frame(y = rnorm(20), x = gl(2, 10))

  ### Ansari-Bradley test
  at <- ansari_test(y ~ x, data = df, distribution = "exact")
  at
  pvalue(at)

  ### bivariate 2-sample problem
  df <- data.frame(y1 = rnorm(20) + rep(0:1, each = 10),
                   y2 = rnorm(20),
                   x = gl(2, 10))

  ### bivariate Pitman-test
  it <- independence_test(y1 + y2 ~ x, data = df,
                          distribution = approximate(B = 9999))
  it
  ### global p-value
  pvalue(it)
  ### joint distribution single-step p-values
  pvalue(it, method = "single-step")
  ### joint distribution (free) step-down p-values
  pvalue(it, method = "step-down")
  ### Bonferroni single-step p-values
  pvalue(it, method = "single-step", distribution = "marginal")
  ### Sidak single-step p-values
  pvalue(it, method = "single-step", distribution = "marginal", type = "Sidak")
  ### Bonferroni (free) step-down p-values
  pvalue(it, method = "step-down", distribution = "marginal")
  ### Sidak (free) step-down p-values
  pvalue(it, method = "step-down", distribution = "marginal", type = "Sidak")
  ### unadjusted p-values
  pvalue(it, method = "unadjusted")

  ### Length of YOY Gizzard Shad from Kokosing Lake, Ohio,
  ### sampled in Summer 1984, Hollander & Wolfe (1999), Table 6.3, page 200
  YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
                               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
                               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
                               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                    site = gl(4, 10, labels = c("I", "II", "III", "IV")))

  ### joint distribution (restricted) step-down p-values for multiple comparisons
  it <- independence_test(length ~ site, data = YOY,
                          xtrafo = mcp_trafo(site = "Tukey"),
                          distribution = approximate(B = 10000))
  it
  pvalue(it, method = "npmcp")

}
\keyword{methods}
\keyword{htest}
