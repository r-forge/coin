% File src/library/stats/man/qqnorm.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2012 R Core Team
% Distributed under GPL 2 or later

\name{qqnorm}
\title{Quantile-Quantile Plots}
\usage{
qqnorm(y, \dots)
\method{qqnorm}{default}(y, ylim, main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE, \dots)

qqline(y, datax = FALSE, distribution = qnorm,
       probs = c(0.25, 0.75), qtype = 7, \dots)

qqplot(x, y, plot.it = TRUE,
       xlab = deparse1(substitute(x)),
       ylab = deparse1(substitute(y)), \dots,
       conf.level = NULL, exact = NULL, 
       simulate = FALSE, B = 2000, col = NA)
}
\alias{qqnorm}
\alias{qqnorm.default}
\alias{qqplot}
\alias{qqline}
\concept{normal probability plot} % PR#12095

\arguments{
  \item{x}{The first sample for \code{qqplot}.}
  \item{y}{The second or only data sample.}
  \item{xlab, ylab, main}{plot labels.  The \code{xlab} and \code{ylab}
    refer to the y and x axes respectively if \code{datax = TRUE}.}
  \item{plot.it}{logical. Should the result be plotted?}
  \item{datax}{logical. Should data values be on the x-axis?}
  \item{distribution}{quantile function for reference theoretical distribution.}
  \item{probs}{numeric vector of length two, representing probabilities.
    Corresponding quantile pairs define the line drawn.}
  \item{qtype}{the \code{type} of quantile computation used in \code{\link{quantile}}.}
  \item{conf.level}{confidence level of the band. The default, \code{NULL},
    does not lead to the computation of a confidence band.}
 \item{exact}{\code{NULL} or a logical indicating whether an exact
    p-value should be computed.  See \sQuote{Details} for the meaning of
    \code{NULL}.}
  \item{simulate.p.value}{a logical indicating whether to compute
    p-values by Monte Carlo simulation.}
  \item{B}{an integer specifying the number of replicates used in the
    Monte Carlo test.}
  \item{col}{the color for filling the confidence band. The default, \code{NA}, is to
          leave the band unfilled.}

  \item{ylim, \dots}{graphical parameters.}
}
\description{
  \code{qqnorm} is a generic function the default method of which
  produces a normal QQ plot of the values in \code{y}.
  \code{qqline} adds a line to a \dQuote{theoretical}, by default
  normal, quantile-quantile plot which passes through the \code{probs}
  quantiles, by default the first and third quartiles.

  \code{qqplot} produces a QQ plot of two datasets. If \code{conf.level} is
  given, a confidence band for a function transforming the distribution of
  \code{x} into the distribution of \code{y} is plotted based on Switzer
  (1976). The QQ plot can be understood as an estimate of such a 
  treatment function. If \code{exact = NULL} (the default), an
  exact confidence band is computed if the product of
  the sample sizes is less than 10000, with or
  without ties. Otherwise, asymptotic distributions are used whose approximations may
  be inaccurate in small samples. Monte-Carlo approximations based on
  \code{B} random permutations are computed when \code{simulate = TRUE}.
  Confidence bands are in agreement with Smirnov's test, that is, the
  bisecting line is covered by the band iff the null of both samples coming 
  from the same distribution cannot be rejected at the same level.

  Graphical parameters may be given as arguments to \code{qqnorm},
  \code{qqplot} and \code{qqline}.
}
\value{
  For \code{qqnorm} and \code{qqplot}, a list with components
  \item{x}{The x coordinates of the points that were/would be plotted}
  \item{y}{The original \code{y} vector, i.e., the corresponding y
    coordinates \emph{including \code{\link{NA}}s}. If \code{conf.level}
  was specified to \code{qqplot}, the list contains additional components 
  \code{lwr} and \code{upr} defining the confidence band.}
}
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988).
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole.

  Switzer, P. (1976).
  Confidence procedures for two-sample problems.
  \emph{Biometrika}, \bold{63}(1), 13--25.
  \doi{10.1093/biomet/63.1.13}.

}
\seealso{
  \code{\link{ppoints}}, used by \code{qqnorm} to generate
  approximations to expected order statistics for a normal distribution.
}
\examples{
require(graphics)

y <- rt(200, df = 5)
qqnorm(y); qqline(y, col = 2)
qqplot(y, rt(300, df = 5))

qqnorm(precip, ylab = "Precipitation [in/yr] for 70 US cities")

## "QQ-Chisquare" : --------------------------
y <- rchisq(500, df = 3)
## Q-Q plot for Chi^2 data against true theoretical distribution:
qqplot(qchisq(ppoints(500), df = 3), y,
       main = expression("Q-Q plot for" ~~ {chi^2}[nu == 3]))
qqline(y, distribution = function(p) qchisq(p, df = 3),
       probs = c(0.1, 0.6), col = 2)
mtext("qqline(*, dist = qchisq(., df=3), prob = c(0.1, 0.6))")
## (Note that the above uses ppoints() with a = 1/2, giving the
## probability points for quantile type 5: so theoretically, using
## qqline(qtype = 5) might be preferable.) 

## Figure 1 in Switzer (1976), knee angle data
switzer <- data.frame(
    angle = c(-31, -30, -25, -25, -23, -23, -22, -20, -20, -18,
              -18, -18, -16, -15, -15, -14, -13, -11, -10, - 9,
              - 8, - 7, - 7, - 7, - 6, - 6, - 4, - 4, - 3, - 2,
              - 2, - 1,   1,   1,   4,   5,  11,  12,  16,  34,
              -31, -20, -18, -16, -16, -16, -15, -14, -14, -14,
              -14, -13, -13, -11, -11, -10, - 9, - 9, - 8, - 7,
              - 7, - 6, - 6,  -5, - 5, - 5, - 4, - 2, - 2, - 2,
                0,   0,   1,   1,   2,   4,   5,   5,   6,  17),
    sex = gl(2, 40, labels = c("Female", "Male")))

ks.test(angle ~ sex, data = switzer)
d <- with(switzer, split(angle, sex))
with(d, qqplot(Female, Male, conf.level = 0.945, col = "lightgrey", 
               exact = TRUE, pch = 19, xlim = c(-31, 31), ylim = c(-31, 31)))
abline(a = 0, b = 1)

## agreement with ks.test
set.seed(1)
x <- rnorm(50)
y <- rnorm(50, mean = .5, sd = .95)
ex <- TRUE
### p = 0.112
(pval <- ks.test(x, y, exact = ex)$p.value)
## 88.8\% confidence band with bisecting line
## touching the lower bound
qqplot(x, y, conf.level = 1 - pval, exact = ex, col = "lightgrey", pch = 19)
abline(a = 0, b = 1)

}
\keyword{hplot}
\keyword{distribution}
