\name{pvalue-methods}
\docType{methods}
\alias{pvalue}
\alias{pvalue-methods}
\alias{pvalue,NullDistribution-method}
\alias{pvalue,IndependenceTest-method}
\alias{pvalue,ScalarIndependenceTest-method}
\alias{pvalue,MaxTypeIndependenceTest-method}
\alias{pvalue,QuadTypeIndependenceTest-method}
\title{ Extract P-Values }
\description{

  Extracts the p-value from objects representing null distributions of 
  independence tests.

}
\usage{
pvalue(object, ...)
}
\arguments{
  \item{object}{an object inheriting from class
    \code{\link{IndependenceTest-class}}.}
  \item{\dots}{additional arguments: \code{method}, 
      a character specifying the type of 
      adjustment (\code{global}, \code{single-step}, \code{step-down}, 
      \code{discrete}, or \code{npmcp})
      should be used. The default is \code{global}.}
}
\section{Methods}{
  \describe{
    \item{pvalue}{extracts the p-value from the specified object.}
  }
}
\details{
  
  Univariate p-values for maximum-type statistics come with associated 99\%
  confidence interval when resampling was used to determine the null
  distribution (which may be the case even when \code{distribution =
  "asypmtotic"} was used). 

  By default, a global p-value is returned. When \code{method =
  "single-step"}, adjusted p-values are obtained from a 
  single-step max-T procedure
  (Westfall & Young, 1993, algorithm 2.5 and formula 2.8). Note that the
  minimum of the adjusted p-values always controls the familywise error 
  rate (FWER) but the maximum type I error, i.e. the error for
  each of the individual tests, is only controlled when the subset
  pivotality condition holds.

  When \code{method = "step-down"} the free step-down resampling method
  (algorithm 2.8 and formula 2.8 in Westfall & Young, 1993) is used, the above
  comments apply as well.

  With \code{method = "discrete"}, the Bonferroni adjustment as suggested by
  Westfall & Wolfinger (1997) with improvements for highly discrete
  permutation distributions is available, however, without taking
  correlations between the test statistics into account. Here, the p-values are
  valid even without assuming subset pivotality.

  For multiple comparisons of levels of one factor (defined via \code{\link{mcp_trafo}}
  p-values can be computed via \code{method = "npmcp"} that are less affected by
  subset pivotality by an algorithm adopted from Westfall (1997).
  This feature is highly experimental.

}
\references{

    Peter H. Westfall \& S. Stanley Young (1993).
    \emph{Resampling-based Multiple Testing}.
    New York: John Wiley & Sons.

    Peter H. Westfall \& Russell D. Wolfinger (1997).
    Multiple tests with discrete distributions. 
    \emph{The American Statistician} \bold{51}, 3--8.

    Peter H. Westfall (1997). Multiple testing of general 
    contrasts using logical constraints and correlation. 
    \emph{Journal of the American Statistical Association}
    \bold{92}, 299--306.

}
\examples{

  ### artificial 2-sample problem
  df <- data.frame(y = rnorm(20), x = gl(2, 10))
 
  ### Ansari-Bradley test
  at <- ansari_test(y ~ x, data = df, distribution = "exact")
  at

  pvalue(at)

  ### bivariate 2-sample problem
  df <- data.frame(y1 = rnorm(20) + c(rep(0, 10), rep(1, 10)), 
                   y2 = rnorm(20), 
                   x = gl(2, 10))

  it <- independence_test(y1 + y2 ~ x, data = df, 
                          distribution = approximate(B = 9999))
  pvalue(it, method = "single-step")
  pvalue(it, method = "step-down")
  pvalue(it, method = "discrete")

  ### Length of YOY Gizzard Shad from Kokosing Lake, Ohio,
  ### sampled in Summer 1984, Hollander & Wolfe (1999), Table 6.3, page 200
  YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
                               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
                               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
                               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                    site = factor(c(rep("I", 10), rep("II", 10),
                                    rep("III", 10), rep("IV", 10))))

  ### p-value for multiple comparisons under logical constraints
  if (require("multcomp")) {
      it <- independence_test(length ~ site, data = YOY,
                              xtrafo = mcp_trafo(site = "Tukey"),
                              distribution = approximate(10000),
                              alternative = "two.sided")
      print(pvalue(it, method = "npmcp"))
  }

}
\keyword{methods}
\keyword{htest}
